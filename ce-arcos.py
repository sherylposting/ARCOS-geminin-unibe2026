# %%
# modules to import
import pandas as pd

# ARCOS
from arcos4py import ARCOS
from arcos4py.tools import filterCollev, calculate_statistics, estimate_eps
from arcos4py.plotting import NoodlePlot, plotOriginalDetrended
from arcos4py import validation
import argparse

# import constants from config.py
from config import *

# %%
# get args for sample output directory from cli
parser = argparse.ArgumentParser()
parser.add_argument('--well', required=True, help='Well ID (e.g., B2)')
parser.add_argument('--site', required=True, help='Site number')
parser.add_argument('--tiffpath', required=True, help='Filepath to the folder of TIFFs')
args = parser.parse_args()

pick_well = args.well
pick_site = args.site
base_dir = args.tiffpath

SAMPLE_OUT_DIR = f"{DIR_OUT_DATA}/well_{pick_well}_site_{pick_site}"

# %%
# Read single-cell time series if saved before
try:
    df_tracks = pd.read_csv(f"{SAMPLE_OUT_DIR}/df_tracks.csv.gz", compression="gzip")
    if not df_tracks.empty:
        print("Successfully loaded previous tracks from df_tracks.csv.gz!")
except Exception as e:
    print(f"Failed to load: {e}")

# Filter out short tracks (shorter than 5 hours)
df_tracks_long = df_tracks[
    df_tracks.groupby(COL_TRK_ID)[COL_TRK_ID].transform("size") >= TRK_MIN_LEN
]

# Create an ARCOS object
obj_arcos = ARCOS(
    data=df_tracks_long,
    position_columns=COL_POS,
    frame_column=COL_FRAME,
    obj_id_column=COL_TRK_ID,
    measurement_column=COL_MEAS,
    clid_column=COL_COLL_ID,
)

# Interpolate and clip time-series data if necessary
if INTERPOLATE:
    obj_arcos.interpolate_measurements()  # interpolation of missing values, optional
if CLIP:
    obj_arcos.clip_meas(CLIP_LOW, CLIP_HIGH)  # clipping of measurement value, optional

# %%
# Detrend and binarize data
obj_arcos.bin_measurements(
    smooth_k=SMOOTHK,
    bias_k=BIASK,
    peak_threshold=PEAK_THRESHOLD,
    binarization_threshold=BINARIZATION_THRESHOLD,
    bias_method=DETRENDING_METHOD,
)  # detrending and binarization

obj_arcos.data.head()

# %%
# detection of CEs
p_nn = estimate_eps(
    obj_arcos.data,
    method="mean",
    position_columns = COL_POS,
    frame_column=COL_FRAME,
    n_neighbors=5,
    plot=True,
    plt_size=(5, 5),
    max_samples=50000,
    mean_multiplier=1
)

# %%
# Detect and track collective events
df_ce = obj_arcos.track_collective_events(
    eps=EPS_COLL, min_clustersize=MIN_CLUSTER_SIZE, n_prev=NPREV_FRAMES
)

# Filter collective events
filterer = filterCollev(
    data=df_ce,
    frame_column=COL_FRAME,
    clid_column=COL_COLL_ID,
    obj_id_column=COL_OBJ_ID,
)  # initialize class

df_ce_filt = filterer.filter(
    min_duration=MIN_COLLECTIVE_EVENT_DURATION, min_total_size=MIN_TOTAL_EVENT_SIZE
)  #

df_ce_filt.head()

# %%
# Calculate statistics of filtered CEs
df_ce_stats = calculate_statistics(
    data=df_ce_filt,
    frame_column=COL_FRAME,
    clid_column=COL_COLL_ID,
    obj_id_column=COL_TRK_ID,
)

median_duration_CE = df_ce_stats["duration"].median()
median_total_size_CE = df_ce_stats["total_size"].median()

mean_duration_CE = df_ce_stats["duration"].mean()
mean_total_size_CE = df_ce_stats["total_size"].mean()

print(df_ce_stats)

if DATA_SAVE:
    df_ce_filt.to_csv(
        f"{SAMPLE_OUT_DIR}/df_ce_filt.csv.gz", compression="gzip", index=False
    )
    obj_arcos.data.to_csv(f"{SAMPLE_OUT_DIR}/df_arcos_data.csv.gz", index=False)  # detrended and binarized, for the time series and eps plots
    df_ce.to_csv(f"{SAMPLE_OUT_DIR}/df_ce.csv.gz", index=False)                   # unfiltered, for the convex hull plot
    df_ce_stats.to_csv(f"{SAMPLE_OUT_DIR}/df_ce_stats.csv.gz", index=False)       # for the size vs duration plot
    print("Successfully saved df_ce_filt.csv.gz, df_arcos_data.csv.gz, df_ce.csv.gz, and df_ce_stats.csv.gz!")