#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 10:15:01 2020

@author: lasse
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import interp1d
import os

import simfun
import sys
import h5py

CCDload = np.load('CCD3D.npy')
QE = np.loadtxt('QE.txt')
optics = np.loadtxt('optics.txt')
spec = np.loadtxt('spec2.txt', skiprows=3)

CCD_size=1000 #subpix
wl_start = 350
wl_stop = 1100
delta_lambda = 1
low_lim=300
up_lim=1150

spectrum = simfun.BMB_interp(x=spec[:,0], y=spec[:,1], wl_start=wl_start, wl_stop=wl_stop, delta_lambda=delta_lambda, kind='cubic', lowlim=low_lim, uplim=up_lim)

QE_int = simfun.BMB_interp(QE[:,0], QE[:,1], wl_start=wl_start, wl_stop=wl_stop, delta_lambda=delta_lambda, kind='cubic', lowlim=low_lim, uplim=up_lim)
optics_int = simfun.BMB_interp(optics[:,0], optics[:,1], wl_start=wl_start, wl_stop=wl_stop, delta_lambda=delta_lambda, kind='quadratic', lowlim=low_lim, uplim=up_lim)
TEC = QE_int[:,1]*optics_int[:,1]
# del QE, optics, spec


#### Generate new PSF: #########
# simfun.psf_maker(res=101, wl_start=350, wl_stop=1100, delta_lambda=1) #can also specify path and filename
# psf_file=simfun.psf_maker(res=101, wl_start=350, wl_stop=1100, delta_lambda=1, si=[120, 120]) #
# psf = psf_file['psf']
# plt.imshow(psf[:,:,0])
################################
file_name = "test_psf2_raw"
file_name = "psf_highpres"
file_name = "psf_ny"
# psf_file = h5py.File(file_name+'.hdf5', 'a') #load in psf_file
# psf = psf_file['psf'] #fetches the psf

psf = np.load(file_name+"_raw.npy")

### Create new 3D CCD? ###
#size = 100
#sub = 10
#full_CCD = np.zeros((size*sub, size*sub, 4))
#for i in np.arange(0,4):
#    full_CCD[:,:,i] = CCD_maker((size, size), subpix=sub) #Will create 3 different CCDs and stack them in a 3D array, 
#    #x,y being the pixel, and z being the "new" ccds
#    sys.stdout.write('*'); sys.stdout.flush(); #"Progress bar", just for visuals 
#np.save('CCD3D.npy', full_CCD)
###
    

ux = int(np.floor(psf.shape[0]/2))
ox = int(np.floor(psf.shape[0]/2)+1)
uy = int(np.floor(psf.shape[1]/2))
oy = int(np.floor(psf.shape[1]/2)+1)
#x_j, y_j = simfun.jitter(steps=1000, gain=0.9, amplitude_act=0.9, amplitude_sens=0.9)
#x_j, y_j = simfun.jitter(steps=1000, gain=0.01, amplitude_act=1.5, amplitude_sens=1.5)
#x_j, y_j = simfun.jitter(steps=10000, gain=1.9, amplitude_act=0.9, amplitude_sens=0.9)
x_j, y_j = simfun.jitter(steps=100, gain=0.2, amplitude_act=3, amplitude_sens=3)


exposure = 100 #Exposure time

###Test of random stars in image###
'''
inp = input('Do you wish to generate random stars? (y/n) ')
if inp == 'y':
    billede=np.zeros((CCD_size,CCD_size))
    #### Exposure time #####
    
    print(' ')
    print('Number of stars complete:')
    print('         10        20        30        40        50        60        70')
    for i in range(30):
        x_pos = np.random.randint(100, 900) 
        y_pos = np.random.randint(100, 900) #generates star position 
        # print(x_pos, ', ', y_pos)
        x_0 = x_pos
        y_0 = y_pos
        
        # plt.plot(x_pos, y_pos, 'r.')
        magni = simfun.mag(np.random.uniform(-1, 6.5)) #Random stellar brightness
        for i in range(exposure):
            billede[x_pos-ux:x_pos+ox, y_pos-uy:y_pos+oy] = billede[x_pos-ux:x_pos+ox, y_pos-uy:y_pos+oy]+psf[:,:,0]*magni #adds psf values to selected area of image array
            x_pos = x_0+x_j[i]
            y_pos = y_0+y_j[i] #updates coordinates based on jitter
            x_pos = int(np.around(x_pos))
            y_pos = int(np.around(y_pos)) # rounds off the coordinates, as matrix can only take int as index
        sys.stdout.write('*'); sys.stdout.flush(); #"Progress bar", just for visuals 
    # image_file = h5py.File('/home/lasse/Documents/uni/Thesis/Sim/image_array.hdf5', "a")
    # image_file.create_dataset('image', shape=(billede.shape[0], billede.shape[1]), dtype='f') #creates datasets in the file.
    # image_file['image'] = billede
        
        if os.path.exists('/home/lasse/Documents/uni/Thesis/Sim/image_file.hdf5') == True: #If file already exists, it will be deleted
            os.remove('/home/lasse/Documents/uni/Thesis/Sim/image_file.hdf5')
        file = h5py.File("image_file.hdf5", "a") #creates the file in which the arrays are stored
        file.create_dataset('image', (CCD_size, CCD_size), dtype='f') #dataset for a CCD image, with stars etc.
        image = file['image']   
        image[:,:] = billede
else: image_file = h5py.File('/home/lasse/Documents/uni/Thesis/Sim/image_file.hdf5', 'a'); billede=image_file['image']
# del ux, ox, uy, oy, x_j, y_j, x_0, y_0,  x_pos, y_pos

#### Slit function

slitwidth= 10
slitheight=150
slitpos = [150, 499]
# slitpos=[499, 499]
mask = simfun.slit(slit_size = [slitwidth, slitheight], pos=slitpos)
billede_masked = billede[:,:]*mask

if os.path.exists("test.fits"):
    os.remove('test.fits')

import astropy.io.fits as fits
hdu = fits.PrimaryHDU(billede_masked)
hdulist = fits.HDUList([hdu])
hdulist.writeto('test.fits')
hdr = hdulist[0].header
# hdr.set('CCD-size', CCD_size)
hdr['CCD-size'] = (str(CCD_size)+'x'+str(CCD_size), 'Dimensions of the CCD detector')
hdr['exptime'] = (exposure, 'Exposure duration (units)')
hdr['slit-w'] = (slitwidth, 'Width of slit')
hdr['slit-h'] = (slitheight, 'Height of slit')
hdr['Slit-x'] = (slitpos[0], 'x-position of slit')
hdr['Slit-y'] = (slitpos[1], 'y-position of slit')


import datetime
day = datetime.date.today()
time = datetime.datetime.now()
hdr.set('date-obs', day.strftime("%m/%d/%y          "))
hdr.set('time-obs', time.strftime("%H:%M:%          "))

hdulist.close()
### Test of random star generation ends here ####
'''
eff = TEC*spectrum[:,1]
slitpos = [150, 499]
billede=np.zeros((CCD_size,CCD_size))
magni = simfun.mag(np.random.uniform(-1, 6.5)) #Random stellar brightness

if os.path.exists("spectrum.fits"): #Iff old files exist, replace by the new one
    os.remove('spectrum.fits')

# hdu = fits.PrimaryHDU(im_disp)
# hdulist = fits.HDUList([hdu])
# hdulist.writeto('spectrum.fits')
# hdr = hdulist[0].header

# # Add information in header
# hdr['CCD-size'] = ('1000'+'x'+str(1000), 'Dimensions of the CCD detector')
# hdr['exptime'] = (exposure, 'Exposure duration (units)')
# hdr['slit-w'] = (slitwidth, 'Width of slit')
# hdr['slit-h'] = (slitheight, 'Height of slit')
# hdr['Slit-x'] = (slitpos[0], 'x-position of slit')
# hdr['Slit-y'] = (slitpos[1], 'y-position of slit')
# hdr['mask'] = ('yes', 'Include slit, y/n')
# hdr['Disp'] = ('yes', 'Include dispersion, y/n')


img_size=[1000, 1000]
# x_dispersion = -1*np.arange(wl_start, wl_stop, delta_lambda)+800
# x_dispersion = np.linspace(0,750, 750)+0.7
# y_dispersion = 0*np.arange(0, 750)
# y_dispersion = 0*np.arange(0, 750)

x_dispersion = 1*np.arange(wl_start, wl_stop, delta_lambda)-400 
# x_dispersion = np.flip(x_dispersion) # flips the x-dispersion
y_dispersion = 0.0002*np.arange(0, 750)**2 # Slight curve in the y-direction
# y_dispersion = 0*np.arange(0, 750)
disper = (x_dispersion, y_dispersion)  #Combine dispersions

# psf = np.flip(psf, axis=2)
eff = np.ones((750))
# mask = simfun.slit(slit_size = [slitwidth, slitheight], pos=slitpos, img_size=[10240, 10240])
mask = simfun.slit(slit_size = [140, 250], pos=slitpos, image_size=img_size) #Create slit-mask
mask=np.ones((img_size)) #Full mask, everything goes through
jit = simfun.jitter_im(x= x_j, y= y_j, psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]) ) #Creating jitter "image"
# jit = np.zeros((101, 101))
# jit[50,50]=1

# psf=np.ones((101, 101))
test, test_wl=simfun.disperser2_copy(wl_endpoints=[350, 1100], jit_img=jit, psf_img=psf, pos=slitpos, image_size=img_size, 
                       dispersion=disper, eff=eff, magni=magni, mask_img=mask, steps=1, plot='n')

ro = simfun.read_out(test)
lam = (test_wl+725*10**(-30))/(test+10**(-30))
# lam = (test_wl)/(test+10**(-30))
plt.figure(); plt.imshow(lam); plt.colorbar()
plt.figure(); plt.plot(ro/ro.max())
# plt.plot(np.arange(92, 750+92), eff/eff.max())




''' for i
i=51

image_size=[1000, 1000]
x_pos=slitpos[0]
y_pos=slitpos[1]
im_disp=np.zeros((1000, 1000))
# psf=np.ones((101,101))
psf=np.zeros((101,101))
psf[50,50] = 1
im_disp=np.zeros((1000, 1000))

for i in range(0,750, steps):
# for i in range(0,101, steps):
    im = np.zeros((image_size[0],image_size[1]))
    fold = simfun.folding(psf[:,:,i], jitter)
    foo = int(psf.shape[0]/2)
    im[0+y_pos-foo:len(fold)+y_pos-foo, 0+x_pos-foo:len(fold)+x_pos-foo] = fold #im[0+y_pos-foo:len(fold)+y_pos-foo, 0+x_pos-foo:len(fold)+x_pos-foo] + fold#*magni
    immask = im*mask
    
    roll_x = np.roll(immask, int(np.modf(  dispersion[i])[1]), axis=1)
    roll_y = np.roll(roll_x,  int(np.modf(y_dispersion[i])[1]), axis=0)

    dx = abs(np.modf(dispersion[i])[0])
    dy = abs(np.modf(y_dispersion[i])[0])
    
    im_disp = im_disp + roll_y*(eff[i]*(1-dx)*(1-dy))  # Add the rolled image to the final, and multiply by the "effectivity"

    roll_dx = np.roll(roll_y, 1, axis=1) # Roll the residual to the next subpixel
    eff_dx = eff[i] * dx * (1-dy) # effectivity of the x-residual
    
    roll_dy = np.roll(roll_y, 1, axis=0) # Roll the residual to the next subpixel, y-wise
    eff_dy = eff[i] * dy * (1-dx) # y-residual eff.
    
    roll_dxy = np.roll(roll_dx, 1, axis=0) # roll the image one step in both x- and y-wise.
    eff_dxy = eff[i]* dx * dy   #and eff.
    
    im_disp = im_disp + roll_dx*eff_dx + roll_dy*eff_dy + roll_dxy*eff_dxy #add all residuals and multiply by their respective effectivities.
            
    sys.stdout.write('/'); sys.stdout.flush(); #"Progress bar", just for visuals     
'''



