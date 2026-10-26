# from config import *


DIR_OUT_DATA = "output-data"
DIR_OUT_PLOT = "output-plots"


IMG_NAPARI = False  # Whether to show images in napari
IMG_NB = True  # Whether to show images in the notebook
DATA_SAVE = True  # Whether to save intermediate data
PLOT_SAVE = True # Whether to save plots
PARCOMP = True  # Whether to perform computations in parallel

if PARCOMP:
    # Parallel processing
    from joblib import Parallel, delayed

    # Number of parallel jobs; set to -1 to use all available cores
    NJOBS = -1


# Analyze only part of the stack (here, chosen frames 10-30)
STACK_TRIM = True

if STACK_TRIM:
    STACK_TRIM_START = 10
    STACK_TRIM_END = 30

# 20-colour tableau palette without greys: #7f7f7f and #c7c7c7
TAB20 = [
    "#1f77b4",
    "#aec7e8",
    "#ff7f0e",
    "#ffbb78",
    "#2ca02c",
    "#98df8a",
    "#d62728",
    "#ff9896",
    "#9467bd",
    "#c5b0d5",
    "#8c564b",
    "#c49c94",
    "#e377c2",
    "#f7b6d2",
    "#bcbd22",
    "#dbdb8d",
    "#17becf",
    "#9edae5",
]

# Minimum nuclei size in pixel (for filtering out small nuclei)
MIN_NUC_SIZE = 200

# Column names in single-cell timeseries from segmentation
COL_POS: list = ["x", "y"]
COL_FRAME: str = "t"
COL_RT: str = "rt"
COL_OBJ_ID: str = "obj_id"
COL_TRK_ID: str = "track_id"
COL_MEAS: str = "cnr"

# column name with collective event track ids
COL_COLL_ID: str = "cl_id"

# Tracking
## Path to btrack config file
BT_CONFIG_FILE: str = "extras/cell_config.json"

## Tracking radius [px]
EPS_TRACK: int = 70

## Minimum track length [frames] (for filtering and passing on to ARCOS)
if STACK_TRIM:
    TRK_MIN_LEN: int = 10
else:
    TRK_MIN_LEN: int = 60

# Time-series pre-processing
## Values for measurement clipping and interpolation
ACQ_FREQ: int = 5
INTERPOLATE: bool = False
CLIP: bool = False
CLIP_LOW: float = 0.001
CLIP_HIGH: float = 0.999

## Binarization parameters
DETRENDING_METHOD: str = "runmed"  # can be "none", "runmed" or "lm"
SMOOTHK: int = 3
BIASK: int = 25  # needed only if BINARIZATION_METHOD method = runmed
PEAK_THRESHOLD: float = 0.1  # needed only if BINARIZATION_METHOD method = runmed
DEGREE_OF_POLYNOMIAL: int = 1  # needed only if BINARIZATION_METHOD method = lm
BINARIZATION_THRESHOLD: float = 0.4

## Collective event detection
EPS_COLL: float = 100
MIN_CLUSTER_SIZE: int = 3
NPREV_FRAMES: int = 1

## Filtering of collective events
MIN_COLLECTIVE_EVENT_DURATION: int = 2
MIN_TOTAL_EVENT_SIZE: int = 3

## Validation of CE detection
N_VALID_STEPS = 200

# crops of a FOV for plotting a single frame
PXSZ = 332.80 / 1024  # pixel size in um

if STACK_TRIM:
    n_frame_show = STACK_TRIM_START # selected frame number to show
else:
    n_frame_show = 47

fov_crop_y = slice(384, 896)  # [px]
fov_crop_x = slice(128, 640)  # [px]