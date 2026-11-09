# %%
import os
import argparse
import numpy as np
import pandas as pd
from skimage.util import map_array

from config import *   # DIR_OUT_DATA, COL_FRAME, COL_OBJ_ID, COL_TRK_ID, COL_COLL_ID, COL_POS, COL_MEAS, TAB20

# %%
parser = argparse.ArgumentParser()
parser.add_argument('--well', required=True, help='Well ID (e.g., B2)')
parser.add_argument('--site', required=True, help='Site number')
parser.add_argument('--tiffpath', required=True, help='Filepath to the folder of TIFFs')
args = parser.parse_args()

SAMPLE_OUT_DIR = f"{DIR_OUT_DATA}/well_{args.well}_site_{args.site}"

# %%
# Load pre-existing data
df_tracks = pd.read_csv(f"{SAMPLE_OUT_DIR}/df_tracks.csv.gz", compression="gzip")
df_ce_filt = pd.read_csv(f"{SAMPLE_OUT_DIR}/df_ce_filt.csv.gz", compression="gzip")
np_lab_filt = np.load(f"{SAMPLE_OUT_DIR}/np_lab_large.npz")["arr"]
np_lab_ring = np.load(f"{SAMPLE_OUT_DIR}/np_lab_ring.npz")["arr"]
n_frames = np_lab_filt.shape[0]
print("Loaded tracks, filtered CEs, nuclear masks and ring masks.")

# %%
# Re-map object IDs to track IDs in labeled images
# Drop untracked objects so NaN track IDs don't end up in the integer label arrays
# (they become background = 0)
df_trk_ok = df_tracks.dropna(subset=[COL_TRK_ID])

np_lab_large_remap = np.zeros_like(np_lab_filt)
np_lab_ring_remap = np.zeros_like(np_lab_ring)
np_lab_large_meas = np.zeros_like(np_lab_filt, dtype="float64")   # cyto/nuc ratio
np_lab_ring_meas = np.zeros_like(np_lab_ring, dtype="float64")

for idx_frame in range(n_frames):
    np_frame = df_trk_ok[df_trk_ok[COL_FRAME] == idx_frame][
        [COL_OBJ_ID, COL_TRK_ID, COL_MEAS]
    ].to_numpy()
    if len(np_frame) == 0:
        continue
    np_nuc_frame = np_lab_filt[idx_frame]
    np_ring_frame = np_lab_ring[idx_frame]
    np_lab_large_remap[idx_frame] = map_array(np_nuc_frame, np_frame[:, 0], np_frame[:, 1])
    np_lab_ring_remap[idx_frame] = map_array(np_ring_frame, np_frame[:, 0], np_frame[:, 1])
    np_lab_large_meas[idx_frame] = map_array(np_nuc_frame, np_frame[:, 0], np_frame[:, 2])
    np_lab_ring_meas[idx_frame] = map_array(np_ring_frame, np_frame[:, 0], np_frame[:, 2])

np.savez_compressed(f"{SAMPLE_OUT_DIR}/np_lab_large_remap.npz", arr=np_lab_large_remap)
np.savez_compressed(f"{SAMPLE_OUT_DIR}/np_lab_ring_remap.npz", arr=np_lab_ring_remap)
np.savez_compressed(f"{SAMPLE_OUT_DIR}/np_lab_large_meas.npz", arr=np_lab_large_meas)
np.savez_compressed(f"{SAMPLE_OUT_DIR}/np_lab_ring_meas.npz", arr=np_lab_ring_meas)
print("Saved track-remapped labels (*_remap.npz) and CNR label images (*_meas.npz)!")

# %%
# CE points
colors = np.take(
    np.array(TAB20),
    [i + 1 for i in np.unique(df_ce_filt[COL_COLL_ID].factorize()[0])],
    mode="wrap",
)

df_w_colors = pd.merge(
    df_ce_filt,
    pd.DataFrame(data={"colors": colors, COL_COLL_ID: df_ce_filt[COL_COLL_ID].unique()}),
    on=COL_COLL_ID,
)

coord_cols = ["axis-0", "axis-1", "axis-2"]   # frame, y, x
out = df_w_colors.rename(columns={
    COL_FRAME: "axis-0",
    COL_POS[1]: "axis-1",
    COL_POS[0]: "axis-2",
})

# hex -> separate r, g, b columns (a tuple in a CSV cell is read back as text)
from matplotlib.colors import to_rgb
out[["r", "g", "b"]] = [to_rgb(c) for c in out["colors"]]

out = out[coord_cols + ["r", "g", "b", COL_COLL_ID]]
out.to_csv(f"{SAMPLE_OUT_DIR}/ce_points.csv", index=False)
print("Saved ce_points.csv!")