# *JWST* MIRI Imaging reduction

Built on [Crab reduction toolkit](https://github.com/1054/Crab.Toolkit.JWST/)

for *JWST* program 3224

with additional intermediate steps by Leonid Sajkov

leonid.sajkov@tufts.edu

Tufts University, June 2024

**Notes**:
- Using the mirim_reduction notebook requires downloading and installing the Crab reduction toolkit (linked above, or at: https://github.com/1054/Crab.Toolkit.JWST/). The notebook was built on toolkit version v20230109_before_huge_mosaic, downloaded May 22, 2024.
- Besides the Crab redution toolkit, this notebook also requries two pieces of code:
    - med_filt_and_sup_bkg_sub.py = code for median filering and super background subtraction;
    - post_reduction_steps.py = code for sorting reduced data and creating TA box cutouts.
- This notebook runs code in the terminal, using the subprocess run module.
