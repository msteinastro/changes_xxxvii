from astropy.io import fits
from astropy.visualization import wcsaxes, simple_norm
import astropy.wcs as wcs
from astropy.wcs import WCS
import astropy.units as u
from regions import Regions
import numpy as np
from scipy.ndimage import label
import matplotlib.pyplot as plt
import plot_constants as pc
import pandas as pd

from photutils.aperture import ApertureMask
from photutils.aperture import aperture_photometry
from photutils.aperture import EllipticalAperture

def edit_header(header):
    header['CRVAL1'] = 0
    header['CRVAL2'] = 0
    return header


def compute_rot_ang(pos_ang, fine_tuning=0):
    print(f"PA={pos_ang} (North Eastward), Fine Tuning={fine_tuning}")
    if pos_ang > 90.:
        print("PA>90: Rotation= -(PA-90)")
        rot_ang =  -(pos_ang - 90.)
        
    else:
        print("PA<90: Rotation= 90-PA")
        rot_ang =  90 - pos_ang
    rot_ang = rot_ang + fine_tuning
    print(f"Rotating by {rot_ang} degrees.")
    # scipy.ndimage.rotate uses an inverse sense of rotation
    rot_ang = -rot_ang
    return rot_ang

def compute_phys_scale(pix_scale, phys_dist):
    """_summary_

    Args:
        pix_scale (_type_): _description_
        phys_dist (_type_): _description_
    """
    phys_pix = np.tan(np.radians(pix_scale))*phys_dist
    phys_pix_kpc = phys_pix.to(u.kpc)
    return phys_pix_kpc

def label_mask(mask):
    fits_mask = fits.open(mask)
    fits_mask.info()
    header = fits_mask[0].header
    mask = np.squeeze(fits_mask[0].data)
    np.savetxt('dummy_mask.csv', mask.astype(int), fmt='%i') # need to save because of big endian error
    new_mask = np.loadtxt('dummy_mask.csv')    
    labeled_mask, n_island = label(new_mask)
    print(f'{n_island} were labeled')
    return labeled_mask   

def my_plot(data, std, mean, wcs_cel, beam, plot_out, exclude=False, title="galaxy"):
    pix_scale = wcs.utils.proj_plane_pixel_scales(wcs_cel)
    beam_artist = beam.ellipse_to_plot(0.1*data.shape[0], 0.1*data.shape[1], np.mean(np.abs(pix_scale)*u.deg))
    beam_artist.set_fill(False)
    beam_artist.set_edgecolor(pc.intensity_overlay)
    beam_artist.set_linestyle('-')
    beam_artist.set_linewidth(0.5)
    norm = simple_norm(data, 'sqrt', min_cut=(mean-5*std), max_cut=0.95*np.nanmax(data))
    fig = plt.figure(figsize=pc.one_c)
    ax = plt.subplot(projection=wcs_cel)
    im = ax.imshow(data, origin='lower', norm=norm, cmap=pc.intensity_cm)
    levels = [3*std *2 ** i for i in range(5)]
    ax.contour(data, levels=levels, colors=pc.intensity_overlay, linewidths=0.5)
    ax.set_ylabel('Declination')
    ax.set_xlabel('Right Ascension')
    ax.set_title(title)
    ax.add_artist(beam_artist)
    if exclude:
        plt.text(0.5, 0.65, 'Excluded', color=pc.intensity_overlay, horizontalalignment='center',
                    verticalalignment='center', transform=ax.transAxes)
    wcsaxes.add_scalebar(ax,1*u.arcmin, label='1arcmin')
    cbar = plt.colorbar(im)
    cbar.set_label('Intensity [Jy/Beam]', horizontalalignment='center')

    plt.axhline(y=data.shape[1]/2, color=pc.intensity_overlay, linestyle='--', linewidth=0.5)
    plt.savefig(plot_out, bbox_inches='tight', dpi=300)
    plt.close()
    
def image_to_dataframe(image_array, pixel_unit, disp_df=True):
    height, width = image_array.shape
    center_x, center_y = width / 2, height / 2
    pixel_values = []
    x_values = []
    z_values = [] 
    for y in range(height):
        for x in range(width):
            if not np.isnan(image_array[y, x]):
                # Calculate Cartesian coordinates relative to center
                delta_x, delta_y = x - center_x, y - center_y
                # Calculate coordinates in the new coordinate system
                x_coord = delta_x 
                z_coord = delta_y 
                rm = image_array[y, x]
                pixel_values.append(rm)
                x_values.append(x_coord)
                z_values.append(z_coord)
    
    df = pd.DataFrame({
        pixel_unit : pixel_values,
        'x': x_values,
        'z': z_values
    })
    # Remove nan's from the dataframe:
    df.dropna(axis=0, inplace=True)
    if disp_df:
        display(df)
    return df

def mask_stacked_cube(rm_map, pi_map,
                      pi_bkg_mean, pi_bkg_std,
                      outname, 
                      use_elliptical_mask=True,
                      nsig_clip=5.0,
                      phy_stack=False):
    pi = fits.getdata(pi_map)
    clip_limit = pi_bkg_mean + nsig_clip*pi_bkg_std
    print(f"Clip all pixels in the RM map where PI<{clip_limit} [Jy/Beam] ({nsig_clip}sigma above mean).")
    pi_mask = np.where(pi>clip_limit, 1. , np.nan)
    rm_fits = fits.open(rm_map)
    rm_head = rm_fits[0].header
    rm_head['COMMENT']="Masked with PI Map (mean+3sig)"
    rm_dat = rm_fits[0].data
    rm_data_masked = rm_dat*pi_mask
    wcs_rm = WCS(rm_head)
    if use_elliptical_mask:
        if phy_stack:
            rm_map_shape = rm_dat.shape
            ellip_ap = EllipticalAperture((rm_map_shape[0]/2, rm_map_shape[1]/2),
                                          a=51,b=25, theta=0)
            reg_mask = ellip_ap.to_mask()
        else:
            reg = Regions.read("rm_elliptical_masks/ang_stacks.reg", format='ds9')
            reg_pix = reg[0].to_pixel(wcs_rm)
            reg_mask = reg_pix.to_mask(mode="center")
        mask_data = reg_mask.to_image(rm_dat.shape)
        mask_data = np.where(mask_data==0, np.nan, 1.0)
        rm_data_masked = rm_data_masked * mask_data
    fits.writeto(outname, data=rm_data_masked, header=rm_head, overwrite=True)
    return 0