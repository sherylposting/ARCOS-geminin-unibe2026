# %% [markdown]
# ## Setup (copied from arcos-tutorial)

# %%
# modules to import
import os
import re
import numpy as np
import pandas as pd
from tifffile import imread, imwrite
from stardist.models import StarDist2D
from csbdeep.utils import normalize
from skimage.morphology import remove_small_objects
import argparse
from glob import glob

# import constants from config.py
from config import *

# %%
# Create directories for outputting plots and data

if not os.path.exists(DIR_OUT_DATA):
    os.makedirs(DIR_OUT_DATA)
    print(f"Directory {DIR_OUT_DATA} created successfully!")
else:
    print(f"Directory {DIR_OUT_DATA} already exists!")

if not os.path.exists(DIR_OUT_PLOT):
    os.makedirs(DIR_OUT_PLOT)
    print(f"Directory {DIR_OUT_PLOT} created successfully!")
else:
    print(f"Directory {DIR_OUT_PLOT} already exists!")

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
# make metadata-labeled output dir
SAMPLE_OUT_DIR = f"{DIR_OUT_DATA}/well_{pick_well}_site_{pick_site}"
os.makedirs(SAMPLE_OUT_DIR, exist_ok=True)

# %%
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

# %% [markdown]
# # Segmentation
# 
# Using a pre-trained [StarDist](https://stardist.net) model `2D_versatile_fluo` for nuclear segmentation.

# %%
# Image stack parameters
n_frames = np_nuc.shape[0]
axis_norm = (0, 1)

print(f"Loaded: {pick_well} Site {pick_site} | {n_frames} frames")

# Create a pretrained model
model = StarDist2D.from_pretrained("2D_versatile_fluo")

# Pre-allocate np array for nuclear labels
np_lab = np.zeros_like(np_nuc)

# Loop over frames, normalize and segment
# Warning: longer execution! Takes approx. 2'30" with GPU support (M2 Max)
for idx_frame in range(n_frames):
    np_frame = np.asarray(np_nuc[idx_frame])
    np_frame = normalize(np_frame, 1.0, 99.8, axis=axis_norm)
    #    labels, details = model.predict_instances(np_frame, prob_thresh=0.5, nms_thresh=0.3)
    ## does it matter to adjust the parameters? why?
    labels, details = model.predict_instances(np_frame)
    np_lab[idx_frame] = labels

# %%
# Filter out small nuclei
# Pre-allocate np array for labels
np_lab_filt = np.zeros_like(np_lab)

for idx_frame in range(n_frames):
    np_lab_filt[idx_frame] = remove_small_objects(
        np_lab[idx_frame],
        min_size=MIN_NUC_SIZE,
    )

# %%
if DATA_SAVE:
    np.savez_compressed(file=f"{SAMPLE_OUT_DIR}/np_lab_large.npz", arr=np_lab_filt)
    print("Successfully saved nuclear masks as np_lab_large.npz!")


