# %%
# modules to import
import os
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from tifffile import imread, imwrite
from skimage.segmentation import expand_labels
from skimage.measure import regionprops_table
import btrack
from btrack.constants import BayesianUpdates
from arcos4py import ARCOS
from arcos4py.tools import filterCollev, calculate_statistics, estimate_eps
import argparse
from glob import glob

# import constants from config.py
from config import *

# %%
# custom functions
# Calculate properties in a labelled region
def loc_calc_regionprops(frame, labels, frameid):
    loc_res = pd.DataFrame(
        regionprops_table(
            labels, frame, properties=["intensity_mean", "centroid", "label"]
        )
    )
    loc_res["frame"] = frameid
    return loc_res

# %%
# get args for data to load from cli
parser = argparse.ArgumentParser()
parser.add_argument('--well', required=True, help='Well ID (e.g., B2)')
parser.add_argument('--site', required=True, help='Site number')
parser.add_argument('--tiffpath', required=True, help='Filepath to the folder of TIFFs')
args = parser.parse_args()

pick_well = args.well
pick_site = args.site
base_dir = args.tiffpath

# %%
SAMPLE_OUT_DIR = f"{DIR_OUT_DATA}/well_{pick_well}_site_{pick_site}"

# Loading image stacks with nuclear and ERK channels
try:
    # Construct filenames from well/site pattern
    h2b_path = os.path.join(base_dir, f"Well{pick_well}_Seq*_Site{pick_site}_C1*.tif")
    gem_path = os.path.join(base_dir, f"Well{pick_well}_Seq*_Site{pick_site}_C2*.tif")
    erk_path = os.path.join(base_dir, f"Well{pick_well}_Seq*_Site{pick_site}_C3*.tif")

    # Use glob if filenames vary
    h2b_path = glob(h2b_path)[0]
    gem_path = glob(gem_path)[0]
    erk_path = glob(erk_path)[0]

    np_nuc = imread(h2b_path)
    np_gem = imread(gem_path)
    np_erk = imread(erk_path)

    if np_nuc.size > 0:
        print(f"Successfully loaded image stacks from {SAMPLE_OUT_DIR}!")
except Exception as e:
    print(f"Failed to load: {e}")

# %%
# If saved before, load labelled nuclear mask from the file
try:
    np_lab_filt = np.load(f"{SAMPLE_OUT_DIR}/np_lab_large.npz", )['arr']
    if np_lab_filt.size > 0:
        print("Successfully loaded previous nuclear masks from np_lab_large.npz!")
except Exception as e:
    print(f"Failed to load: {e}")

n_frames = np_lab_filt.shape[0]

if PARCOMP:
    # Process labels in parallel
    def process_frame(idx_frame):
        return expand_labels(np_lab_filt[idx_frame], 6) - expand_labels(np_lab_filt[idx_frame], 2) 

    np_lab_ring = np.array(
        Parallel(n_jobs=NJOBS)(delayed(process_frame)(idx_frame) for idx_frame in range(n_frames))
    )
else:
    # Process labels sequentially by frame
    np_lab_ring = np.zeros_like(np_lab_filt)

    # Loop over frames to filter out small nuclei
    for idx_frame in range(n_frames):
        np_lab_ring[idx_frame] = expand_labels(np_lab_filt[idx_frame], 6) - expand_labels(np_lab_filt[idx_frame], 2)

# Extract measurements from images based on nuclear and ring masks.
# Extracted: intensity_mean, centroid, label
PARCOMP = True
if PARCOMP:
    # Calculate properties in parallel
    df_prop_nuc = pd.concat(
        Parallel(n_jobs=NJOBS)(
            delayed(loc_calc_regionprops)(img, lab, fid)
            for (fid, (img, lab)) in enumerate(zip(np_erk, np_lab_filt))
        )
    )
    df_prop_cyto = pd.concat(
        Parallel(n_jobs=NJOBS)(
            delayed(loc_calc_regionprops)(img, lab, fid)
            for (fid, (img, lab)) in enumerate(zip(np_erk, np_lab_ring))
        )
    )

else:
    # Calculate properties sequentially by frame
    df_prop_nuc = pd.DataFrame()
    df_prop_cyto = pd.DataFrame()

    for idx_frame in range(n_frames):
        df_prop_nuc_frame = pd.DataFrame(
            regionprops_table(
                np_lab_filt[idx_frame],
                np_erk[idx_frame],
                properties=["intensity_mean", "centroid", "label"],
            )
        )

        df_tmp_prop_cyto = pd.DataFrame(
            regionprops_table(
                np_lab_ring[idx_frame], np_erk[idx_frame], properties=["intensity_mean", "label"]
            )
        )
        df_prop_nuc_frame["frame"] = idx_frame
        df_tmp_prop_cyto["frame"] = idx_frame

        df_prop_nuc = pd.concat([df_prop_nuc, df_prop_nuc_frame])
        df_prop_cyto = pd.concat([df_prop_cyto, df_tmp_prop_cyto])

# Rename columns names and merge into a single df with properties
df_prop_nuc.rename(
    columns={
        "frame": COL_FRAME,
        "label": COL_OBJ_ID,
        "intensity_mean": "intensity_mean_nuc",
        "centroid-0": COL_POS[1],
        "centroid-1": COL_POS[0],
    },
    inplace=True,
)

df_prop_cyto.rename(
    columns={
        "frame": COL_FRAME,
        "label": COL_OBJ_ID,
        "intensity_mean": "intensity_mean_cyto",
        "centroid-0": COL_POS[1],
        "centroid-1": COL_POS[0],
    },
    inplace=True,
)

# Check whether columns exist before dropping
if df_prop_cyto.columns.isin(COL_POS).any():
    df_prop_cyto.drop(columns=COL_POS, inplace=True)

#df_prop = df_prop_nuc.merge(df_prop_cyto, how="inner", on=[COL_FRAME, COL_OBJ_ID])
df_prop = pd.merge(df_prop_nuc, 
                   df_prop_cyto, 
                   how="inner", 
                   on=[COL_FRAME, COL_OBJ_ID])

df_prop.head()

# Create btrack objects (with properties) from the segmentation data
# By setting properties=('area',) the tracking algorithm takes into account the mask's area in addition to x/y.
# By setting assign_class_ID=True adds the object ID that will be used later in merging tracking results with the df_prop table.
obj_bt = btrack.utils.segmentation_to_objects(
    np_lab_filt, properties=("area",), assign_class_ID=True
)

# Initialise a tracker session using a context manager
with btrack.BayesianTracker() as tracker:
    # configure the tracker using a config JSON file
    tracker.configure(BT_CONFIG_FILE)

    tracker.update_method = BayesianUpdates.APPROXIMATE
    tracker.max_search_radius = EPS_TRACK

    # append the objects to be tracked
    tracker.append(obj_bt)

    # set the volume (Z axis volume limits default to [-1e5, 1e5] for 2D data)
    tracker.volume = ((0, np_nuc.shape[2]), (0, np_nuc.shape[1]))

    # track them (in interactive mode)
    tracker.track_interactive(step_size=100)

    # generate hypotheses and run the global optimizer
    tracker.optimize()

    # get the tracks as a python list
    ls_tracks = tracker.tracks

    # optional: get the data in a format for napari
    bt_data, bt_props, bt_graph = tracker.to_napari()

# Convert btrack output to pandas
df_tracks_bt = pd.concat(
    pd.DataFrame(t.to_dict(["ID", "t", "x", "y"])) for t in ls_tracks
)

# Rename columns to make them consistent with the rest
df_tracks_bt.rename(
    columns={
        "ID": COL_TRK_ID,
        "t": COL_FRAME,
        "x": COL_POS[0],
        "y": COL_POS[1],
        "class_id": COL_OBJ_ID,
    },
    inplace=True,
)

# Set column types
df_tracks_bt[COL_OBJ_ID] = df_tracks_bt[COL_OBJ_ID].astype(
    "Int32"
)  # Int32/64 allows integer nan's https://stackoverflow.com/a/65580754/1898713

# Merge btrack output with the initial single-cell data from segmentation
df_tracks = pd.merge(
    df_prop[
        [
            COL_FRAME,
            COL_OBJ_ID,
            COL_POS[0],
            COL_POS[1],
            "intensity_mean_nuc",
            "intensity_mean_cyto",
        ]
    ],
    df_tracks_bt[[COL_FRAME, COL_OBJ_ID, COL_TRK_ID]].dropna(),
    on=[COL_FRAME, COL_OBJ_ID],
    how="left",
)

# Calculate cyto/nuc ratio
df_tracks[COL_MEAS] = df_tracks["intensity_mean_cyto"] / df_tracks["intensity_mean_nuc"]

# Add a column with real time in hours
df_tracks[COL_RT] = df_tracks[COL_FRAME] * ACQ_FREQ / 60

# final df with track id's
df_tracks.head()

if DATA_SAVE:
    df_tracks.to_csv(f"{SAMPLE_OUT_DIR}/df_tracks.csv.gz", compression="gzip", index=False)
    np.savez_compressed(f"{SAMPLE_OUT_DIR}/np_lab_ring.npz", arr=np_lab_ring)
    print(f"Successfully saved tracks (with properties) as df_tracks.csv.gz, and labeled cyto rings as np_lab_ring.npz!")