#### Input parameters ####

in_spec = '../sample_values/star1_vega_nm.dat' 		#Input science spectrum
in_spec2 = '../sample_values/star1_vega_nm_w_exoplanet.dat' 		#Input science spectrum

col_area = 200. 			# Collecting area

img_size = [100, 1024] 			#Size of the CCD, in pixels

pl_scale = 348.3664772727273	# Plate scale "/mm

pix_size = 0.0135			# mm/pixel 

bg_spec= '../sample_values/zodiac.txt' 		# Background spectrum, i.e. zodiacal light

exp = 300 				# Exposure time

sub_pixel = 10 			# Amount of sub-pixels per full pixel

wl_ran = [300, 1000] 			# Wavelength range

eta_in = '../sample_values/QE2.txt', '../sample_values/optics.txt'		# Spectral troughput of the entire system. Requires at minimum the CCD QE

slit = ['pix', 2, 3.5] 		# Slit size. Unit first, then width and height

#psf = 'psf_stor_raw.npy' 			# Point Spread Function of the optics etc.
# psf = ''

psf_col = ''                    #PSF colours

disper = '../sample_values/disp.npy' 		#Dispersion of the spectrograph

####### Optionals ########
jitter = '' 				# Spacecraft jitter

step = 10 				# Step size. Only needed if jitter is left 'default'

in_CCD = '../sample_values/CCD.npy' 		# Input CCD imperfections. Sub-pixel variations

CCD_col = [300.,  600.,  850., 1100.] # CCD colours, respective to each slice in in_CCD

