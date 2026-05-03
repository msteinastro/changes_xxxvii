import astropy.units as u
tar_res_scale = 115 * u.arcsec
cal_distance_Mpc = 10
regrid_pix_scale_arcsec = 5 # 2.5
regrid_npix = 600 # 1200
n_diam_ir = 5
max_ang_extent_arcmin = 11.53
excluded_galaxies = ['NGC660','NGC4438', 'NGC4594', 'NGC5084', 'UGC10288', # Exculuded in Intensity stack (Wiegert+2015)
                     'NGC2992'
                    ]
# angular stacking
regrid_npix_ang = 150
ang_stack_pix_per_beam = 5
ang_stack_inner_radius_beams = 2.5
ang_stack_outer_radius_beams = 7.5
ang_stack_min_sampling = 15.0
# highres angular stacking
regrid_npix_ang_hr = 250            
ang_stack_pix_per_beam_hr = 7.0
ang_stack_min_sampling_hr = 19.0

high_res_sample_ang = ["NGC4631","NGC4565","NGC3628",
                       "NGC891","NGC5907","NGC4192","NGC3556",
                       "NGC2683","NGC2613","NGC3079","NGC4666",
                       "NGC4096","NGC4217","NGC4302","NGC5775"]
# physical stacking
phy_res_kpc = 2.0
phy_pix_size_kpc = phy_res_kpc/5
extent_phy_ref_frame_kpc = 50
n_pix_phy_ref_frame = int(extent_phy_ref_frame_kpc/phy_pix_size_kpc)
ref_distance_Mpc = 20
phy_stack_inner_radius_kpc = 5
phy_stack_outer_radius_kpc = 12
# high res stacking
phy_pix_size_kpc_hr = 0.2
extent_phy_ref_frame_kpc_hr = 50
n_pix_phy_ref_frame_hr = int(extent_phy_ref_frame_kpc_hr/phy_pix_size_kpc_hr)
# highres sample physical resolution
high_res_sample = ['NGC2683','NGC4631','NGC3628','NGC891', 'NGC3432',
                   'NGC4096','NGC4565','NGC4192', 'NGC3556','NGC4157',
                   'NGC4013','NGC4388','NGC5907', 'NGC3877']
#stacking rm_maps of individual galaxies:
min_n_for_stack=10
abs_rm_val_norm = 50
# added systematic RM uncertatinty derived from synthetic data:
sys_rm_err = 90 # rad/m^2