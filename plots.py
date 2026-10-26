#!/usr/bin/env python3
"""Plotting code from arcos-tutorial.ipynb, for one well/site. Loads saved results only."""

import os
import argparse
from glob import glob

import numpy as np
import pandas as pd
from tifffile import imread

import matplotlib
matplotlib.use("Agg")  # no display on cluster nodes
from matplotlib import pyplot as plt

from scipy.spatial import ConvexHull
import plotnine as p9
from mizani.formatters import percent_format

from skimage.util import montage

from arcos4py.tools import estimate_eps
from arcos4py.plotting import NoodlePlot, plotOriginalDetrended

from config import *

# %%
# args
parser = argparse.ArgumentParser()
parser.add_argument("--well", required=True, help="Well ID (e.g., B2)")
parser.add_argument("--site", required=True, help="Site number")
parser.add_argument("--tiffpath", type=str, required=True, help="Folder of TIFFs")
args = parser.parse_args()

pick_well = args.well
pick_site = args.site
base_dir = args.tiffpath

SAMPLE_OUT_DIR = f"{DIR_OUT_DATA}/well_{pick_well}_site_{pick_site}"
DIR_OUT_PLOT = f"{DIR_OUT_PLOT}/well_{pick_well}_site_{pick_site}"
os.makedirs(DIR_OUT_PLOT, exist_ok=True)

# %%
try:
    # Load raw images
    np_nuc = imread(glob(os.path.join(base_dir, f"Well{pick_well}_Seq*_Site{pick_site}_C1*.tif"))[0])
    np_erk = imread(glob(os.path.join(base_dir, f"Well{pick_well}_Seq*_Site{pick_site}_C3*.tif"))[0])

    # Load saved results
    np_lab_filt = np.load(f"{SAMPLE_OUT_DIR}/np_lab_large.npz")["arr"]
    np_lab_ring = np.load(f"{SAMPLE_OUT_DIR}/np_lab_ring.npz")["arr"]

    df_tracks = pd.read_csv(f"{SAMPLE_OUT_DIR}/df_tracks.csv.gz", compression="gzip")
    df_arcos_data = pd.read_csv(f"{SAMPLE_OUT_DIR}/df_arcos_data.csv.gz", compression="gzip")
    df_ce = pd.read_csv(f"{SAMPLE_OUT_DIR}/df_ce.csv.gz", compression="gzip")
    df_ce_filt = pd.read_csv(f"{SAMPLE_OUT_DIR}/df_ce_filt.csv.gz", compression="gzip")
    df_ce_stats = pd.read_csv(f"{SAMPLE_OUT_DIR}/df_ce_stats.csv.gz", compression="gzip")

    print(f"Successfully loaded image stacks and outputs for {SAMPLE_OUT_DIR}!")
except Exception as e:
    print(f"Failed to load: {e}")

# Trim the stack
if STACK_TRIM:
    np_nuc = np_nuc[35:55]
    np_erk = np_erk[35:55]

# %%
# Show a montage of 3 time points from the nuclear and ERK channels
if IMG_NB:
    fig, axarr = plt.subplots(2, 1)

    np_nuc_sel = np_nuc[np.linspace(0, len(np_nuc) - 1, 3).astype(np.int32)]
    np_erk_sel = np_erk[np.linspace(0, len(np_erk) - 1, 3).astype(np.int32)]

    imMontNuc = montage(
        np_nuc_sel,
        padding_width=10,
        fill=np.max(np_nuc),
        grid_shape=(1, np_nuc_sel.shape[0]),
    )

    imMontErk = montage(
        np_erk_sel,
        padding_width=10,
        fill=np.max(np_erk),
        grid_shape=(1, np_erk_sel.shape[0]),
    )

    axarr[0].set_axis_off()
    axarr[0].set_title("Nuclear channel: H2B")
    axarr[0].imshow(imMontNuc, cmap="gray", vmin=50, vmax=8000)

    axarr[1].set_axis_off()
    axarr[1].set_title("ERK channel: ERK-KTR")
    axarr[1].imshow(imMontErk, cmap="gray", vmin=500, vmax=7000)

    if PLOT_SAVE:
        plt.savefig(f"{DIR_OUT_PLOT}/montage-raw.png", bbox_inches="tight", dpi=300)
    plt.close()

# %%
np_nuc_show = np_nuc[n_frame_show, fov_crop_y, fov_crop_x]

plt.imshow(np_nuc_show, cmap="Greys", interpolation="none", vmin=1000, vmax=5000)
plt.axis("off")
plt.tight_layout()

if PLOT_SAVE:
    plt.savefig(
        f"{DIR_OUT_PLOT}/fig2-raw-nuc.png", bbox_inches="tight", pad_inches=0.0, dpi=300
    )
plt.close()

# %%
np_erk_show = np_erk[n_frame_show, fov_crop_y, fov_crop_x]

plt.imshow(np_erk_show, cmap="Greys", interpolation="none", vmin=1500, vmax=7000)
plt.axis("off")
plt.tight_layout()

if PLOT_SAVE:
    plt.savefig(
        f"{DIR_OUT_PLOT}/fig2-raw-erk.png", bbox_inches="tight", pad_inches=0.0, dpi=300
    )
plt.close()

# %%
# Convert zeros to nan's, so the background is displayed in white
np_lab_show = np.where(
    np_lab_filt[n_frame_show, fov_crop_y, fov_crop_x] == 0,
    np.nan,
    np_lab_filt[n_frame_show, fov_crop_y, fov_crop_x],
)

plt.imshow(np_lab_show, cmap="tab20", interpolation="none")
plt.axis("off")
plt.tight_layout()

if PLOT_SAVE:
    plt.savefig(
        f"{DIR_OUT_PLOT}/fig2-seg-nuc.png", bbox_inches="tight", pad_inches=0.0, dpi=300
    )
plt.close()

# %%
# Overlay nuclear masks on EKR-KTR
plt.imshow(np_nuc_show, cmap="Greys", interpolation="none", vmin=1000, vmax=5000)
plt.imshow(np_lab_show, cmap="tab20", alpha=0.25, interpolation="none")
plt.axis("off")
plt.tight_layout()

if PLOT_SAVE:
    plt.savefig(
        f"{DIR_OUT_PLOT}/fig3-raw-nuc-seg.png", bbox_inches="tight", pad_inches=0.0, dpi=300
    )
plt.close()

# %%
np_lab_ring_show = np.where(
    np_lab_ring[n_frame_show, fov_crop_y, fov_crop_x] == 0,
    np.nan,
    np_lab_ring[n_frame_show, fov_crop_y, fov_crop_x],
)

plt.imshow(np_lab_ring_show, cmap="tab20", interpolation="none")
plt.axis("off")
plt.tight_layout()

if PLOT_SAVE:
    plt.savefig(
        f"{DIR_OUT_PLOT}/fig2-seg-ring.png", bbox_inches="tight", pad_inches=0.0, dpi=300
    )
plt.close()

# %%
# Overlay cytoplasmic rings on EKR-KTR
plt.imshow(np_erk_show, cmap="Greys", interpolation="none", vmin=1500, vmax=7000)
plt.imshow(np_lab_ring_show, cmap="tab20", alpha=0.75, interpolation="none")
plt.axis("off")
plt.tight_layout()

if PLOT_SAVE:
    plt.savefig(
        f"{DIR_OUT_PLOT}/fig3-raw-erk-ring.png",
        bbox_inches="tight",
        pad_inches=0.0,
        dpi=300,
    )
plt.close()

# %%
# Both labels
plt.imshow(np_lab_show, cmap="tab20", interpolation="none")
plt.imshow(np_lab_ring_show, cmap="tab20", alpha=0.75, interpolation="none")
plt.axis("off")
plt.tight_layout()

if PLOT_SAVE:
    plt.savefig(
        f"{DIR_OUT_PLOT}/fig3-seg-nuc-ring.png",
        bbox_inches="tight",
        pad_inches=0.0,
        dpi=300,
    )
plt.close()

# %%
# 3D projection of all identified tracks
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")
ax.view_init(120, 150, 61.5)
ax.set_proj_type("persp", focal_length=0.4)
ax.set_box_aspect(aspect=(1, 1, 3))
ax.invert_xaxis()
ax.invert_yaxis()
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_zlabel("")
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

v_track_ids = df_tracks[COL_TRK_ID].unique()
v_colors = plt.cm.tab20(np.linspace(0, 1, len(v_track_ids)))

for track_idx, track_id in enumerate(v_track_ids):
    df_track = df_tracks[df_tracks[COL_TRK_ID] == track_id]
    ax.plot(
        df_track[COL_POS[0]],
        df_track[COL_POS[1]],
        df_track[COL_RT],
        color=v_colors[track_idx],
        linewidth=0.5,
    )

if PLOT_SAVE:
    plt.savefig(
        f"{DIR_OUT_PLOT}/fig2-tracks3D.png", bbox_inches="tight", pad_inches=0.0, dpi=300
    )
plt.close()

# %%
# Random sample of 200 tracks
np.random.seed(3)
v_track_ids_rnd = np.random.choice(df_tracks["track_id"].unique(), 200)

s_tmp_x = f"{COL_FRAME} * {ACQ_FREQ} / 60"
s_tmp_y = f"{COL_POS[1]} * {PXSZ}"

p_tracks_proj = (
    p9.ggplot(
        df_tracks[df_tracks[COL_TRK_ID].isin(v_track_ids_rnd)],
        p9.aes(x=s_tmp_x, y=s_tmp_y, group=COL_TRK_ID),
    )
    + p9.geom_line(p9.aes(color="pd.Categorical(track_id)"), size=0.5, alpha=0.8)
    + p9.scale_color_manual(
        values=np.take(np.array(TAB20), range(0, len(v_track_ids_rnd)), mode="wrap")
    )
    + p9.xlab("Time [hours]")
    + p9.ylab("Position Y [μm]")
    + p9.theme_bw()
    + p9.theme(
        legend_position="none",
        panel_grid=p9.element_blank(),
        panel_border=p9.element_rect(size=0.5, colour="black"),
        panel_background=p9.element_blank(),
        text=p9.element_text(),
        axis_text_x=p9.element_text(colour="black", size=10),
        axis_text_y=p9.element_text(colour="black", size=10),
        axis_title_x=p9.element_text(colour="black", size=12),
        axis_title_y=p9.element_text(colour="black", size=12),
    )
)

if PLOT_SAVE:
    p9.ggsave(p_tracks_proj,
              f"{DIR_OUT_PLOT}/fig5-tracks-proj-y.pdf",
              width=3,
              height=2)

# %%
# Track lengths
df_tracks_len = (
    df_tracks.groupby("track_id")["t"]
    .apply(lambda x: x.max() - x.min() + 1)
    .reset_index(name="track_len")
)
df_tracks_len["track_len"] = df_tracks_len["track_len"] * ACQ_FREQ / 60

p_tracks_len_hist = (
    p9.ggplot(df_tracks_len, p9.aes(x="track_len", y=p9.after_stat("width*density")))
    + p9.geom_histogram(binwidth=1, boundary=0, color="black")
    + p9.geom_vline(
        xintercept=TRK_MIN_LEN * ACQ_FREQ / 60, color="red", linetype="dashed"
    )
    + p9.scale_y_continuous(labels=percent_format())
    + p9.xlab("Track length [hours]")
    + p9.ylab("width x density")
    + p9.theme_bw()
    + p9.theme(
        legend_position="none",
        panel_grid=p9.element_blank(),
        panel_border=p9.element_rect(size=0.5, colour="black"),
        panel_background=p9.element_blank(),
        text=p9.element_text(),
        axis_text_x=p9.element_text(colour="black", size=10),
        axis_text_y=p9.element_text(colour="black", size=10),
        axis_title_x=p9.element_text(colour="black", size=12),
        axis_title_y=p9.element_text(colour="black", size=12),
    )
)

if PLOT_SAVE:
    p9.ggsave(
        p_tracks_len_hist,
        f"{DIR_OUT_PLOT}/fig5-tracks-hist-len.pdf",
        width=2.5,
        height=2
    )

# %%
plt.imshow(np_lab_show, cmap="viridis", interpolation="none")
plt.axis("off")
plt.tight_layout()

if PLOT_SAVE:
    plt.savefig(
        f"{DIR_OUT_PLOT}/fig2-seg-nuc-erk.png", bbox_inches="tight", pad_inches=0.0, dpi=300
    )
plt.close()

# %%
# Sample ERK activity time series (df_arcos_data holds the tracks kept for ARCOS)
np.random.seed(7)
v_track_ids_rnd = np.random.choice(df_arcos_data["track_id"].unique(), 10)

s_tmp_x = f"{COL_FRAME} * {ACQ_FREQ} / 60"
s_tmp_y = f"{COL_MEAS}"

p_tracks_erk = (
    p9.ggplot(
        df_arcos_data[df_arcos_data[COL_TRK_ID].isin(v_track_ids_rnd)],
        p9.aes(x=s_tmp_x, y=s_tmp_y, group=COL_TRK_ID),
    )
    + p9.geom_line(p9.aes(color="pd.Categorical(track_id)"), size=0.5, alpha=0.8)
    + p9.scale_color_manual(
        values=np.take(np.array(TAB20), range(0, len(v_track_ids_rnd)), mode="wrap")
    )
    + p9.xlab("Time [hours]")
    + p9.ylab("Position Y [μm]")
    + p9.theme_bw()
    + p9.theme(
        legend_position="none",
        panel_grid=p9.element_blank(),
        panel_border=p9.element_rect(size=0.5, colour="black"),
        panel_background=p9.element_blank(),
        text=p9.element_text(),
        axis_text_x=p9.element_text(colour="black", size=10),
        axis_text_y=p9.element_text(colour="black", size=10),
        axis_title_x=p9.element_text(colour="black", size=12),
        axis_title_y=p9.element_text(colour="black", size=12),
    )
)

if PLOT_SAVE:
    p9.ggsave(p_tracks_erk, f"{DIR_OUT_PLOT}/fig5-tracks-erk.pdf", width=3, height=2)

# %%
p_detrended = plotOriginalDetrended(
    df_arcos_data,
    frame_column = COL_FRAME,
    measurement_column = COL_MEAS,
    detrended_column = "cnr.resc",
    obj_id_column = COL_TRK_ID,
    seed=13,
)

p_detrended.plot_original_and_detrended(n_samples=25, subplots=(5, 5), plotsize=(20, 10), add_binary_segments = True)
if PLOT_SAVE:
    plt.savefig(f"{DIR_OUT_PLOT}/fig6-tracks-detrended.png", bbox_inches="tight", dpi=300)

fig, ax = p_detrended.plot_original(
    n_samples=9, subplots=(3, 3), plotsize=(9, 5), add_binary_segments=True
)

# Custom labels for x-axis ticks
labels = np.arange(0, 360, step=60) * ACQ_FREQ / 60
labels = labels.astype(int)

for a in ax.flatten():
    a.set_xticks(np.arange(0, 360, step=60))
    a.set_xticklabels(labels, )

fig.supxlabel("Time [hours]")
fig.supylabel("C/N ERK-KTR")

if PLOT_SAVE:
    fig.savefig(f"{DIR_OUT_PLOT}/fig6-tracks-bin.svg", bbox_inches="tight", pad_inches=0.0)
plt.close("all")

# %%
# eps estimation
p_nn = estimate_eps(
    df_arcos_data,
    method="mean",
    position_columns = COL_POS,
    frame_column=COL_FRAME,
    n_neighbors=5,
    plot=True,
    plt_size=(5, 5),
    max_samples=50000,
    mean_multiplier=1
)
if PLOT_SAVE:
    plt.savefig(f"{DIR_OUT_PLOT}/fig6-estimate-eps.png", bbox_inches="tight", dpi=300)
plt.close("all")

# %%
# Noodle plot
fig, ax = NoodlePlot(
    df_ce_filt[df_ce_filt["cl_id"] == 5], COL_COLL_ID, COL_TRK_ID, COL_FRAME, COL_POS[0], COL_POS[1]
).plot("y", color_cylce=TAB20)

v_lab_x = np.arange(0, 360, step=60) * ACQ_FREQ / 60
v_lab_x = v_lab_x.astype(int)

ax.set_xticks(np.arange(0, 360, step=60))
ax.set_xticklabels(v_lab_x, )

v_lab_y = np.arange(0, 1024, step=100 / PXSZ) * PXSZ
v_lab_y = v_lab_y.astype(int)

ax.set_yticks(np.arange(0, 1024, step=100 / PXSZ))
ax.set_yticklabels(v_lab_y, )

ax.invert_yaxis()

plt.xlabel("Time [hours]")
plt.ylabel("Position Y [μm]")

if PLOT_SAVE:
    fig.savefig(f"{DIR_OUT_PLOT}/fig6-noodle.svg", bbox_inches="tight", pad_inches=0.0)
plt.close("all")

# %%
# Plot size vs duration of CEs
median_duration_CE = df_ce_stats["duration"].median()
median_total_size_CE = df_ce_stats["total_size"].median()

mean_duration_CE = df_ce_stats["duration"].mean()
mean_total_size_CE = df_ce_stats["total_size"].mean()

p_ce_dur_sz = (
    p9.ggplot(df_ce_stats, p9.aes(x="total_size", y="duration * ACQ_FREQ"))
    + p9.geom_point(alpha=0.1)
    + p9.geom_point(
        data=pd.DataFrame(
            data={
                "total_size": [median_total_size_CE],
                "duration": [median_duration_CE],
            }
        ),
        color="red",
        size=3,
    )
    + p9.xlab("Total number of cells involved in CE")
    + p9.ylab("Duration of CE [minutes]")
    + p9.theme_bw()
    + p9.theme(
        legend_position="right",
        panel_grid=p9.element_blank(),
        panel_border=p9.element_rect(size=0.5, colour="black"),
        panel_background=p9.element_blank(),
        text=p9.element_text(),
        axis_text_x=p9.element_text(colour="black", size=10),
        axis_text_y=p9.element_text(colour="black", size=10),
        axis_title_x=p9.element_text(colour="black", size=12),
        axis_title_y=p9.element_text(colour="black", size=12),
    )
)

if PLOT_SAVE:
    p9.ggsave(
        p_ce_dur_sz, f"{DIR_OUT_PLOT}/fig6-ce-dur-sz.pdf", width=3.5, height=3, dpi=300
    )

# %%
# Sample CE: convex hull on ERK-KTR image
plt.close()

df_ce_show = df_ce[
    (df_ce[COL_FRAME] == n_frame_show)
    & (df_ce[COL_POS[0]].between(fov_crop_x.start, fov_crop_x.stop))
    & (df_ce[COL_POS[1]].between(fov_crop_y.start, fov_crop_y.stop))
]

m_xy_ce = df_ce_show[[COL_POS[0], COL_POS[1]]].values
chull = ConvexHull(m_xy_ce)

np_erk_show = np_erk[n_frame_show, fov_crop_y, fov_crop_x]
plt.imshow(np_erk_show, cmap="Greys", interpolation="none", vmin=1500, vmax=7000)

for simplex in chull.simplices:
    plt.plot(
        m_xy_ce[simplex, 0] - fov_crop_x.start,
        m_xy_ce[simplex, 1] - fov_crop_y.start,
        "-",
        color="#ff7f0e",
        linewidth=4,
    )
plt.plot(
    m_xy_ce[:, 0] - fov_crop_x.start,
    m_xy_ce[:, 1] - fov_crop_y.start,
    "o",
    color="#d62728",
)

plt.axis("off")
plt.tight_layout()

if PLOT_SAVE:
    plt.savefig(
        f"{DIR_OUT_PLOT}/fig2-ce-erk.png", bbox_inches="tight", pad_inches=0.0, dpi=300
    )
plt.close()