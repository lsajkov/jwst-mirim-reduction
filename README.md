# *JWST* MIRI Imaging reduction

Built on Crab reduction toolkit for JWST program 3224
with additional intermediate steps by Leonid Sajkov (leonid.sajkov@tufts.edu)

Tufts University, June 2024

Notes:
- Using the mirim_reduction.ipynb notebook requires downloading and installing the Crab reduction toolkit (https://github.com/1054/Crab.Toolkit.JWST/). The notebook was built on toolkit version v20230109_before_huge_mosaic, downloaded May 22, 2024.
- Besides the Crab redution toolkit, this notebook also requries two pieces of code:
    - med_filt_and_sup_bkg_sub.py = code for median filering and super background subtraction;
    - post_reduction_steps.py = code for sorting reduced data and creating TA box cutouts.
- This notebook runs code in the terminal, using the subprocess run module.

Pipeline steps:
Where necessary, user input is listed. Data reduction parameters can be calibrated in step 8.

### 1 Import necessary directories

### 2 Set necessary paths/parameters
Inputs:
glbl_dir - directory for all final project files
loc_dir - local directory for storing intermediate steps (files can get very large, so syncing them with cloud servers can take a while)

### SAVE LOGS: Y/N
# IMPORTANT: the text log outputs from this process are long and can take a lot of memory
# when running a notebook like this one locally, caching the logs (as many editors do automatically)
# can cause the editor to crash/slow down the entire process. Change to true to save output logs. 
save_outputs = False 

### 3 Login and download data from MAST
Inputs:
program_id - JWST program ID, for downloading files from MAST
MAST_login_token - MAST login token, for downloading files from MAST

### 4 Open working directory and copy uncal files

### 5 Precache necessary CRDS files

### 6 Stage 1: uncal -> rate

### 7 Stage 2.1: rate -> cal
SKIP_BKGSUB - toggle to skip uniform background subtraction step using skymatch. Additional subtraction will be done after stage 3, pass 1
DETECT_THRESH - equivalent to SExtractor parameter DETECT_THRESH
OVERWRITE - toggle to overwrite prior generated files from this 

### Optional: make backups of cal files:

### 8 Intermediate steps
Inputs:
MEDIAN_FILTER = toggle to perform median filering on the "hot" pixels in the image (default: True)
WINDOW_SIZE = if performing median filtering, what size the "window" should be. Takes int values. (default: 1)
SUPER_BKG_SUB = toggle to perform super background subtraction on the images (default: True)
SMALL_SOURCE_THRESH = remove small sources from the segmentation map, so that they will be considered in the super median background. Use to avoid "hot" pixels and spurious detections clogging the background. If 0, does not remove small sources. Takes int values. (default: 0)
DILATE_SEG = toggle to dilate the segmentation map when masking sources for super background.
DILATE_ITERATIONS = how many iterations of the seg map dilation to perform. Takes int values. (default: 2)
MEDIAN_THRESH = minimum number of images in a FILTER/DURATION stack in order to perform subtraction e.g.: if median-thresh = 5 and there are 3 exposures in F1500W lasting 88.8s, none of these images will have a super median background subtracted. Takes int values. (default: 5)
SUPPRESS_WARNINGS = suppress default python warnings, as functions often encounter all-nan slices that they feel the need to warn us about

### 9 Stage 2.2: cal -> i2d

### 10 Package results/produce TA box
Inputs:
DESTINATION_DIRECTORY = where to place reduced, packaged files
WORKING_DIRECTORY = the reducing_data folder, if not already specified in above cells
DATE = date the reduction was done, for bookkeeping. Optional. Takes string. (default: '', defaults to today's date)
MAKE_CUTOUTS = toggle to make cutouts of the TA box (produces both pdf and fits)

Readme made on March 24, 2025
Contact Leonid Sajkov at leonid.sajkov@tufts.edu with questions, comments, corrections, and bugs.
