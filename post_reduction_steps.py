#!/usr/bin/env python
#
"""""
Post reduction steps in MIRI imaging background reduction

Written by Leonid Sajkov
Tufts University, June 2024
leonid.sajkov@tufts.edu

Input:

Process:

Output:
"""""
import os, click, fnmatch, shutil
from datetime import datetime

import pysiaf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel


def get_header_info(file):
    with fits.open(file) as fitsfile:

        global sci_img, wcs, header
        
        sci_img = fitsfile[1].data
        wcs = WCS(fitsfile[1].header)
        header = fitsfile[1].header

        global RA_V1, DEC_V1, PA_V3, TARGPROP, OBSLABEL, OBSLABEL_concat, FILTER

        RA_V1    = fitsfile[1].header['RA_V1']
        DEC_V1   = fitsfile[1].header['DEC_V1']
        PA_V3    = fitsfile[1].header['PA_V3']
        TARGPROP = fitsfile[0].header['TARGPROP']

        OBSLABEL = fitsfile[0].header['OBSLABEL']
        FILTER = fitsfile[0].header['FILTER']

        OBSLABEL_concat = ''.join(OBSLABEL.split(' '))

def isolate_TA_LRS_ROI():
    
    global sky_ra, sky_dec
    global cntrCoord, cntrCoord_pixel

    siaf = pysiaf.Siaf('MIRI')
    ta = siaf['MIRIM_TALRS']
    print('{} V2Ref = {}'.format(ta.AperName, ta.V2Ref))
    print('{} V3Ref = {}'.format(ta.AperName, ta.V3Ref))
    for attribute in ['InstrName', 'AperShape']:
        print('{} {} = {}'.format(ta.AperName, attribute, getattr(ta, attribute)))

    sci_x = ta.corners('sci')[0]
    sci_y = ta.corners('sci')[1]

    attmat = pysiaf.utils.rotations.attitude_matrix(0, 0, RA_V1, DEC_V1, PA_V3)
    ta.set_attitude_matrix(attmat)
    sky_ra, sky_dec = ta.sci_to_sky(sci_x, sci_y)
    ta_roi_cntr_ra, ta_roi_cntr_dec = ta.tel_to_sky(ta.V2Ref, ta.V3Ref) 

    cntrCoord = SkyCoord(ra = ta_roi_cntr_ra, dec = ta_roi_cntr_dec, unit = 'deg')
    cntrCoord_pixel = skycoord_to_pixel(cntrCoord, wcs)

def make_cutout(D = 60):

    global cutout
    cutout = Cutout2D(sci_img, cntrCoord, 2 * D, wcs = wcs)

def plot_ROI(D = 60,
             vmin = -1, vmax = 4, save = False, savename = ''):

    fig = plt.figure(figsize = (10, 10))
    ax = fig.add_subplot(projection = cutout.wcs)

    ax.imshow(cutout.data, origin = 'lower', vmin = vmin, vmax = vmax)
    ax.scatter(np.mean(sky_ra), np.mean(sky_dec), transform = ax.get_transform('fk5'), s = 550, marker = 'x', c = 'r')


    p = patches.Polygon(list(zip(sky_ra, sky_dec)), edgecolor = 'r',
        facecolor = 'none', transform = ax.get_transform('fk5'), lw = 3)
    
    ax.text(0.01, 0.01, f'{OBSLABEL}    {np.mean(sky_ra):.3f}{np.mean(sky_dec):+.3f}',
             transform = ax.transAxes, color = 'white', fontsize = 20)
        
    ax.add_patch(p)
    ax.axis('off')

    if save:
        if savename != '':
            fig.savefig(savename, dpi = 150,
                         bbox_inches = 'tight')
        else:
            raise(NameError('Provide a destination for the cutout'))

def save_fits_file(fitspath):

    fits_file = fits.ImageHDU()
    fits_file.header = header
    fits_file.header.update(cutout.wcs.to_header())
    fits_file.data = cutout.data
    fits_file.writeto(fitspath, overwrite = True)

@click.command()
@click.argument('destination_directory', type = click.Path())
@click.argument('working_directory', default = '', type = click.Path())
@click.option('--date', default = '')
@click.option('--make-cutouts/--no-make-cutouts', is_flag = True, default = True)

def main(
         destination_directory,
         working_directory,
         date,
         make_cutouts
    ):

    if working_directory == '':
        working_directory = os.getcwd()

    click.echo('Looking for a calibrated3_mosaics subdirectory...')
    if ~working_directory.endswith('calibrated3_mosaics'):
        working_directory += '/calibrated3_mosaics'
    
    click.echo(f'Working directory: {working_directory}')

    if not os.path.exists(working_directory): raise(FileNotFoundError('calibrated3_mosaics directory not found.'))

    all_files = fnmatch.filter(os.listdir(working_directory), 'jw*_obs*_MIRI_*')
    if len(all_files) == 0: raise(FileNotFoundError('No jw*_obs*_MIRI_* directories found in working directory.'))

    if date == '':
        date_today = datetime.today()
        date = date_today.strftime('%d%b%Y')

    for file in all_files:

        with fits.open(f'{working_directory}/{file}/{file}_i2d.fits') as hdu:
            i2d_header = hdu[0].header
            OBSLABEL = i2d_header['OBSLABEL']
            OBSNUMBER = i2d_header['OBSERVTN']
            TARGPROP = i2d_header['TARGPROP']

        FILTER = file.split('_')[-1]
        
        if not os.path.exists(f'{destination_directory}/{TARGPROP}'):
            os.mkdir(f'{destination_directory}/{TARGPROP}')
        
        if not os.path.exists(f'{destination_directory}/{TARGPROP}/{FILTER}'):
            os.mkdir(f'{destination_directory}/{TARGPROP}/{FILTER}')

        shutil.copy(f'{working_directory}/jw03224_obs{OBSNUMBER}_MIRI_{FILTER}/jw03224_obs{OBSNUMBER}_MIRI_{FILTER}_i2d.fits',
                    f'{destination_directory}/{TARGPROP}/{FILTER}/{TARGPROP}_MIRI_{FILTER}_{date}_i2d.fits')
        
        if make_cutouts:

            get_header_info(f'{working_directory}/{file}/{file}_i2d.fits')

            isolate_TA_LRS_ROI()

            make_cutout()

            plot_ROI(save = True,
                     savename = f'{destination_directory}/{TARGPROP}/{FILTER}/{TARGPROP}_MIRI_{FILTER}_{date}_TA_cutout.pdf')
            
            save_fits_file(f'{destination_directory}/{TARGPROP}/{FILTER}/{TARGPROP}_MIRI_{FILTER}_{date}_TA_cutout.fits')

if __name__ == '__main__':
    main()