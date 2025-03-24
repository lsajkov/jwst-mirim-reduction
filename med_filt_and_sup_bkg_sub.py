#!/usr/bin/env python
#
"""""
Intermediate steps in MIRI imaging background reduction

Written by Leonid Sajkov
Tufts University, June 2024
leonid.sajkov@tufts.edu

Input:
- MIRI images reduced up to stage 3 (i2d) that have NOT been rescaled (rescaling pixel scale = 1)
- corresponding segmentation maps
- corresponding stage 2 (cal) files
- files must be of the form: ./calibrated3_mosaics/jw*_obs*_MIRI_*/jw*_obs*_MIRI_*_i2d.fits (or *segm.fits)
  AND the corresponding _cal.fits files, in the same directory (./calibrated3_mosaics/jw*_obs*_MIRI_*/jw*mirimage_cal.fits)
Process:
- Median filters images with a given kernel (size and type specified using kernel-size and kernel-type)
- Creates a super median background from the given dataset:
  - groups together observations of the same filter and same exposure time
  - subtracts them from correspodning filtered images.
Output:
- 

"""""

import os
import click

import cv2
import shutil
import fnmatch

import numpy as np
from astropy.io import fits

# Define functions
def median_filter_image(cal_img, seg_map, window_size):

    #fil nan values with zeroes, to avoid later confusion with masked pixels
    nan_pixel_mask = np.where(np.isnan(cal_img),True, False)
    seg_map_mask = np.where(seg_map == 0, True, False)

    row_len, col_len = np.shape(cal_img)

    filtered_cal_img = cal_img

    ws = window_size

    #median filter by row
    for ii in range(row_len):
        for jj in range(col_len):
            if nan_pixel_mask[ii, jj] & seg_map_mask[ii, jj]:
                slice = cal_img[ii - ws: ii + ws + 1,
                                jj - ws: jj + ws + 1]
                filtered_cal_img[ii, jj] = np.nanmedian(slice)

    return filtered_cal_img

# Main
@click.command()
@click.argument('working_directory', default = '', type = click.Path())
@click.option('--median-filter/--no-median-filter', is_flag = True, default = True)
@click.option('--window-size', type = int, default = 1)
@click.option('--super-bkg-sub/--no-super-bkg-sub', is_flag = True, default = True)
@click.option('--small-source-thresh', type = int, default = 0)
@click.option('--dilate-seg/--no-dilate-seg', is_flag = True, default = True)
@click.option('--dilate-iterations', type = int, default = 2)
@click.option('--median-thresh', type = int, default = 5)
@click.option('--debug-skip-filtering', is_flag = True, default = False)
@click.option('--debug-skip-bkg-sub', is_flag = True, default = False)
def main(
         working_directory,
         median_filter,
         window_size,
         super_bkg_sub,
         small_source_thresh,
         dilate_seg,
         dilate_iterations,
         median_thresh,
         debug_skip_filtering,
         debug_skip_bkg_sub
    ):
    
    if working_directory == '':
        working_directory = os.getcwd()
    
    click.echo(f"Median filtering:\t{'YES' if (median_filter) & (~debug_skip_filtering)  else 'NO'}")    
    click.echo(f"Super background sub:\t{'YES' if (super_bkg_sub) & (~debug_skip_bkg_sub)  else 'NO'}")

    all_files = fnmatch.filter(os.listdir(working_directory), 'jw*_mirimage')
    # median filter all applicable files
    i = 0
    if (median_filter) & (~debug_skip_filtering):
        click.echo('Median filtering...')
        click.echo(f'\r[{int(30 * i/len(all_files)) * "*"}{(30 - int(30 * i/len(all_files))) * "-"}] ({i}/{len(all_files)})', nl = False)
        for file in all_files:

            files_here = os.listdir(f'{working_directory}/{file}/calibrated2_cals')

            cal_file = fnmatch.filter(files_here, '*_cal.fits')

            if len(cal_file) == 1:
                with fits.open(f"{working_directory}/{file}/calibrated2_cals/{cal_file[0]}") as hdu:
                    cal_img = hdu[1].data
            else: raise(NameError('Too many/not enough cal files. There should only be one in the directory.'))

            seg_file = f'{working_directory}/{file}/calibrated2_cals/{file}_cal_run_sextractor_classic_dir/SExtractor_Segmentation.fits'
            if os.path.exists(seg_file):
                with fits.open(seg_file) as hdu:
                    seg_map = hdu[0].data
            else: raise(NameError('Too many/not enough seg files. There should only be one in the directory.'))                   

            filtered_cal_img = median_filter_image(cal_img = cal_img, seg_map = seg_map, window_size = window_size)
                
            #save filtered image as an intermediate '_filt.fits' file
            with fits.open(f"{working_directory}/{file}/calibrated2_cals/{cal_file[0]}") as hdu:

                hdu[1].data = filtered_cal_img

                filt_path  = cal_file[0].split('_')
                filt_path[-1] = 'filt.fits'
                filt_path = '_'.join(filt_path)
                hdu.writeto(f"{working_directory}/{file}/calibrated2_cals/{filt_path}", overwrite = True)

            i += 1
            click.echo(f'\r[{int(30 * i/len(all_files)) * "*"}{(30 - int(30 * i/len(all_files))) * "-"}] ({i}/{len(all_files)})', nl = False)
        
        click.echo('\nMedian filtering complete\n')

    if (super_bkg_sub) & (~debug_skip_bkg_sub):
        click.echo('Super background subtracting...')
        #First, construct superstack of backgrounds:
        #   dictionary of 2-d (x/y) arrays labeled by:
        #       1. FILTER
        #       2. Exposure time (DURATION)
        #       3. Exposure number

        SUPERSTACK = {}

        for file in all_files:

            with fits.open(f"{working_directory}/{file}/calibrated2_cals/{file}_cal.fits") as hdu:
                cal_header = hdu[0].header
                FILTER = cal_header['FILTER']
                DURATION = cal_header['DURATION']
                duration_key = ''.join(str(DURATION).split('.'))[:-2]

            if f"{FILTER}_{duration_key}" not in SUPERSTACK.keys():
                SUPERSTACK[f"{FILTER}_{duration_key}"] = [] #creates "open spot" in the superstack for the image

            if median_filter:

                if os.path.exists(f"{working_directory}/{file}/calibrated2_cals/{file}_filt.fits"):
                    with fits.open(f"{working_directory}/{file}/calibrated2_cals/{file}_filt.fits") as hdu:
                        image = hdu[1].data
                else: raise(NameError('Filt file does not exist.'))

            elif ~median_filter:
                with fits.open(f"{working_directory}/{file}/calibrated2_cals/{file}_cal.fits") as hdu:
                    image = hdu[1].data

            seg_file = f'{working_directory}/{file}/calibrated2_cals/{file}_cal_run_sextractor_classic_dir/SExtractor_Segmentation.fits'
            if os.path.exists(seg_file):
                with fits.open(seg_file) as hdu:
                    seg_map = hdu[0].data
            else: raise(NameError('Too many/not enough seg files. There should only be one in the directory.'))

            #remove smaller emission
            if small_source_thresh > 0:
                for source in range(np.nanmax(seg_map)):
                    if np.sum(seg_map == source) <= small_source_thresh:
                        seg_map = np.where(seg_map == source, 0, seg_map)

            #dilate seg map to account for extended emission:
            if dilate_seg:
                dilation_kernel = np.ones((3, 3), np.uint8)       
                seg_map = np.array(np.where(seg_map == 0, 0, 1), dtype = np.uint8)
                seg_map = cv2.dilate(seg_map, dilation_kernel, iterations = dilate_iterations)

            masked_image = np.where(seg_map == 0, image, np.nan)
            SUPERSTACK[f"{FILTER}_{duration_key}"].append(masked_image/np.nanmedian(masked_image)) #get median to be 1.0

        
        # Make median background for each filter, exposure pair
        MEDIAN_SUPERSTACK = {}
        for key in SUPERSTACK:

            MEDIAN_SUPERSTACK[key] = np.nanmedian(SUPERSTACK[key], axis = 0)

        for key in MEDIAN_SUPERSTACK:
            MEDIAN_SUPERSTACK[key] = MEDIAN_SUPERSTACK[key]/np.nanmedian(MEDIAN_SUPERSTACK[key])
            if (np.nanmedian(MEDIAN_SUPERSTACK[key]) - 1) > 0.05: raise(ValueError(f'Median for {key} not 1.0. It is {np.nanmedian(MEDIAN_SUPERSTACK[key])}'))

        np.save('superstack.py', SUPERSTACK)
        np.save('median_superstack.npy', MEDIAN_SUPERSTACK)

        click.echo('Super median background created.\nSubtracting...')
        # We now have a dictionary labeled by FILTER_EXPTIME that contains a median background for each avaialble filter/exptime pair

        # Now subtract
        i = 0
        click.echo(f'\r[{int(30 * i/len(all_files)) * "*"}{(30 - int(30 * i/len(all_files))) * "-"}] ({i}/{len(all_files)})', nl = False)
        for file in all_files:
            if median_filter:
                
                with fits.open(f"{working_directory}/{file}/calibrated2_cals/{file}_filt.fits") as hdu:
                    filt_header = hdu[0].header
                    FILTER = filt_header['FILTER']
                    DURATION = filt_header['DURATION']
                    duration_key = ''.join(str(DURATION).split('.'))[:-2]

                    #will only do median subtraction if there are at least as many images as a defined threshold (default: 5)
                    if len(SUPERSTACK[f'{FILTER}_{duration_key}']) >= median_thresh:
                        filt_img = hdu[1].data

                        median_bkg = MEDIAN_SUPERSTACK[f'{FILTER}_{duration_key}']
                        median_bkg *= np.nanmedian(filt_img)/np.nanmedian(median_bkg)

                        subbed_img = filt_img - median_bkg
                        hdu[1].data = subbed_img
                        hdu.writeto(f"{working_directory}/{file}/calibrated2_cals/{file}_sub.fits",
                                    overwrite = True)
                    
                    else:
                        hdu.writeto(f"{working_directory}/{file}/calibrated2_cals/{file}_sub.fits",
                                    overwrite = True)


            elif ~median_filter:

                seg_file = f'{working_directory}/{file}/calibrated2_cals/{file}_cal_run_sextractor_classic_dir/SExtractor_Segmentation.fits'
                if os.path.exists(seg_file):
                    with fits.open(seg_file) as hdu:
                        seg_map = hdu[0].data
                else: raise(NameError('Too many/not enough seg files. There should only be one in the directory.'))   

                with fits.open(f"{working_directory}/{file}/calibrated2_cals/{file}_cal.fits") as hdu:
                    cal_header = hdu[0].header
                    FILTER = cal_header['FILTER']
                    DURATION = cal_header['DURATION']
                    duration_key = ''.join(str(DURATION).split('.'))[:-2]

                    #will only do median subtraction if there are at least as many images as a defined threshold (default: 5)
                    if len(SUPERSTACK[f'{FILTER}_{duration_key}']) >= median_thresh: 
                        cal_img = hdu[1].data

                        median_bkg = MEDIAN_SUPERSTACK[f'{FILTER}_{duration_key}']
                        median_bkg *= np.nanmedian(cal_img)

                        subbed_img = cal_img - median_bkg
                        
                        hdu[1].data = subbed_img
                        hdu.writeto(f"{working_directory}/{file}/calibrated2_cals/{file}_sub.fits",
                                    overwrite = True)
                    
                    else:
                        hdu.writeto(f"{working_directory}/{file}/calibrated2_cals/{file}_sub.fits",
                                    overwrite = True)

            i += 1
            click.echo(f'\r[{int(30 * i/len(all_files)) * "*"}{(30 - int(30 * i/len(all_files))) * "-"}] ({i}/{len(all_files)})', nl = False)
        
        click.echo('\nBackground subtraction complete\n')

    for file in all_files:
        
        seg_file = f'{working_directory}/{file}/calibrated2_cals/{file}_cal_run_sextractor_classic_dir/SExtractor_Segmentation.fits'
        if os.path.exists(seg_file):
            with fits.open(seg_file) as hdu:
                seg_map = hdu[0].data
        else: raise(NameError('Too many/not enough seg files. There should only be one in the directory.'))

        if os.path.exists(f"{working_directory}/{file}/calibrated2_cals/{file}_sub.fits"):

            with fits.open(f"{working_directory}/{file}/calibrated2_cals/{file}_sub.fits") as hdu:
                subbed_img = hdu[1].data

            with fits.open(f"{working_directory}/{file}/calibrated2_cals/{file}_cal.fits") as hdu:
                cal_img = hdu[1].data

                hdu[1].data = subbed_img
                hdu.writeto(f"{working_directory}/{file}/calibrated2_cals/{file}_cal.fits", overwrite = True)

        elif os.path.exists(f"{working_directory}/{file}/calibrated2_cals/{file}_filt.fits"):
            click.echo('Saving cal from filt files')
            
            with fits.open(f"{working_directory}/{file}/calibrated2_cals/{file}_filt.fits") as hdu:
                filt_img = hdu[1].data

            with fits.open(f"{working_directory}/{file}/calibrated2_cals/{file}_cal.fits") as hdu:
                cal_img = hdu[1].data

                hdu[1].data = np.where(seg_map == 0, filt_img, cal_img)
                hdu.writeto(f"{working_directory}/{file}/calibrated2_cals/{file}_cal.fits", overwrite = True)

if __name__ == '__main__':
    main()