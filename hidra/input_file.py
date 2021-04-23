#### Input parameters ####

in_spec = '../sample_values/star1_vega_nm.dat' 		#Input science spectrum. Must be 2 columns, first containing wavelengths, second containing the photons/s/cm^2/Å
in_spec2 = '../sample_values/star1_vega_nm_w_exoplanet.dat' 		#Input science spectrum, same format as in_spec

col_area = 200. 			# Collecting area, in cm^2

img_size = [100, 1024] 			#Size of the CCD, in no. of pixels

pl_scale = 348.3664772727273	# Plate scale, in "/mm

pix_size = 0.0135			# pixel size, in mm/pixel 

bg_spec= '../sample_values/zodiac.txt' 		# Background spectrum, i.e. zodiacal light. Same unit as in_spec

exp = 300 				# Exposure time, in seconds

sub_pixel = 10 			# Amount of sub-pixels per full pixel

wl_ran = [300, 1000] 			# Wavelength range, in nm

eta_in = '../sample_values/QE2.txt', '../sample_values/optics.txt'		# Spectral troughput of the entire system. Requires at minimum the CCD QE. Values in the arrays must be between 0 and 1

slit = ['pix', 2, 3.5] 		# Slit size. Unit first, then width and height measured in

#psf = 'psf_stor_raw.npy' 			# Point Spread Function of the optics etc. Is now obsolete. 
# psf = ''

psf_col = ''                    #PSF colours (wavelength endpoints), in nm

disper = '../sample_values/disp.npy' 		#Dispersion of the spectrograph. Two coloumns, first is the dispersion in x, measured in pixels, second column is dispersion in y, also in pixels. There must be one value per wavelength

####### Optionals ########
jitter = '' 				# Spacecraft jitter

step = 10 				# Number of steps per second exposure

in_CCD = '../sample_values/CCD.npy' 		# Input CCD imperfections. Sub-pixel variations. 2D array, with values close to 1. The further from 1, the more "imperfect" a detector.

CCD_col = [300.,  600.,  850., 1100.] # CCD colours, respective to each slice in in_CCD, in nm. Not useful as of the current version 

