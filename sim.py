#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 10:15:01 2020
Hyperchromatic instrument spectra simulator 
  HISS
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

import matplotlib
textwidth=337
height = 3
width = textwidth/72.27
fontsize = 11

plt.style.use('default')
# latex_preamble = [
#     r'\usepackage{lmodern}',
#     r'\usepackage{amsmath}',
#     r'\usepackage{amsfonts}',
#     r'\usepackage{mathtools}',
# ]
matplotlib.rcParams.update({
    # 'text.usetex'        : True,
    # 'font.family'        : 'sans',
#    'font.serif'         : 'cmr10',
    'font.size'          : fontsize,
    # 'mathtext.fontset'   : 'cm'#,
#    'text.latex.preamble': latex_preamble,
})


CCDload = np.load('CCD3D.npy')
wls = np.linspace(350, 1100, 4)

QE = np.loadtxt('QE.txt')
optics = np.loadtxt('optics.txt')
# spec = np.loadtxt('spec2.txt', skiprows=3)
spec = np.loadtxt('spec3.txt', skiprows=4)
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
# file_name = "test_psf2_raw"
# file_name = "psf_highpres"
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
  
'''

ux = int(np.floor(psf.shape[0]/2))
ox = int(np.floor(psf.shape[0]/2)+1)
uy = int(np.floor(psf.shape[1]/2))
oy = int(np.floor(psf.shape[1]/2)+1)
#x_j, y_j = simfun.jitter(steps=1000, gain=0.9, amplitude_act=0.9, amplitude_sens=0.9)
#x_j, y_j = simfun.jitter(steps=1000, gain=0.01, amplitude_act=1.5, amplitude_sens=1.5)
#x_j, y_j = simfun.jitter(steps=10000, gain=1.9, amplitude_act=0.9, amplitude_sens=0.9)



expo = 100 #Exposure time in seconds
exp_step_size = 10 # num. of steps per second

exposure = int(expo*exp_step_size) #exposure, steps

spectrum[:,1] = spectrum[:,1]*(exp_step_size**-1)

x_j, y_j = simfun.jitter(steps=exposure, gain=0.2, amplitude_act=3, amplitude_sens=3)



eff = TEC*spectrum[:,1]
slitpos = [150, 499]
billede=np.zeros((CCD_size,CCD_size))
magni = simfun.mag(np.random.uniform(-1, 6.5)) #Random stellar brightness

# if os.path.exists("spectrum.fits"): #Iff old files exist, replace by the new one
#     os.remove('spectrum.fits')

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


# x_dispersion = np.flip(x_dispersion) # flips the x-dispersion
# y_dispersion = 0*np.arange(0, 750)
x_dispersion = 1*np.arange(wl_start, wl_stop, delta_lambda)-400 
y_dispersion = 0.0002*np.arange(0, 750)**2 # Slight curve in the y-direction
disper = (x_dispersion, y_dispersion)  #Combine dispersions

# psf = np.flip(psf, axis=2)
# eff = np.ones((750))
# mask = simfun.slit(slit_size = [slitwidth, slitheight], pos=slitpos, img_size=[10240, 10240])
mask = simfun.slit(slit_size = [140, 250], pos=slitpos, image_size=img_size) #Create slit-mask
# mask=np.ones((img_size)) #Full mask, everything goes through
jit = simfun.jitter_im(x= x_j, y= y_j, psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]) ) #Creating jitter "image"
# jit = np.zeros((101, 101))
# jit[50,50]=1

# psf=np.ones((101, 101))
test, test_wl=simfun.disperser2(wl_endpoints=[350, 1100], jit_img=jit, psf_img=psf, pos=slitpos, image_size=img_size, 
                       dispersion=disper, eff=eff, magni=magni, mask_img=mask, steps=1, plot='n')

# ro = simfun.read_out(test)
lam = (test_wl+725*10**(-80))/(test+10**(-80))
# lam = (test_wl)/(test+10**(-30))
# plt.figure(); plt.imshow(lam); plt.colorbar() 
np.save("img_wl.npy", lam)
np.save("img.npy", test)
# plt.figure(); plt.plot(ro/ro.max())
# plt.plot(np.arange(92, 750+92), eff/eff.max())

CCDload = np.load('CCD3D.npy')
wls = np.linspace(350, 1100, 4)
img = np.load('img.npy')
img_wl = np.load('img_wl.npy')
t =simfun.ccd_interp(CCDload, wls, test, lam)
ro = simfun.read_out(t)

binned = simfun.binner(t, 10)
ro = simfun.read_out(binned)

# plt.figure(); plt.imshow(t2)

# final = simfun.ccd_interp(CCDload, wls, test, test_wl)


# plt.figure()
# plt.imshow(t)
# plt.xlabel('Sub-pixel')
# plt.ylabel('Sub-pixel')

ticksiz = 14
fontsiz = 15
# plt.figure(figsize=(6, 5))
plt.plot(np.arange(0,750)+98, eff*exp_step_size, label='Input spectrum'); 
plt.plot(ro, label=str(expo)+' second exposure')
# plt.title(str(expo)+' second exposure, with ' + str(exp_step_size)+ ' steps per second')
plt.ylabel('Photons', size=14)
plt.xlabel('Sub-pixel and $\lambda$', size=14)
plt.grid()
# plt.yticks(fontsize= ticksiz); plt.yticks(fontsize= ticksiz)
plt.legend()
'''


img_size=[1000, 1000] # in sub-pixels
expo = 100 #Exposure time in seconds
exp_step_size = 10 # num. of steps per second
exposure = int(expo*exp_step_size) #exposure, steps

spectrum[:,1] = spectrum[:,1]*(exp_step_size**-1)
eff = TEC*spectrum[:,1]

slitpos = [150, 499]

x_dispersion = 1*np.arange(wl_start, wl_stop, delta_lambda)-400 
y_dispersion = 0.0002*np.arange(0, 750)**2 # Slight curve in the y-direction
disper = (x_dispersion, y_dispersion)  #Combine dispersions
mask = simfun.slit(slit_size = [140, 250], pos=slitpos, image_size=img_size)
magni = simfun.mag(np.random.uniform(-1, 6.5))

CCDload = np.load('CCD3D.npy')
wls = np.linspace(350, 1100, 4)

NN = int(input('How many runs? '))
read_outs = np.zeros((100, NN))
for i in range(NN):
    x_j, y_j = simfun.jitter(steps=exposure, gain=0.2, amplitude_act=3, amplitude_sens=3)

    jit = simfun.jitter_im(x= x_j, y= y_j, psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]) ) #Creating jitter "image"
    image, image_wl=simfun.disperser2(wl_endpoints=[350, 1100], jit_img=jit, psf_img=psf, pos=slitpos, image_size=img_size, 
                       dispersion=disper, eff=eff, magni=magni, mask_img=mask, steps=1, plot='n')
    lam = (image_wl+725*10**(-80))/(image+10**(-80))

    image_ccd = simfun.ccd_interp(CCDload, wls, image, lam)
    
    binned = simfun.bin_sum(image_ccd, 10)
    stoej = np.random.normal(0, 0.05, [100,100])*np.sqrt(sum(sum(binned)))
    
    ro = simfun.read_out(binned+stoej)
    read_outs[:,i] = ro
    
sum(eff*10)*100/(sum(sum(binned+stoej))/np.mean(CCDload))
from PIL import Image
import simfun
image = np.array(Image.open('coffee.jpg').convert('L'))
image_bin = simfun.bin_sum(image, 20)
plt.figure(figsize=[4.8, 4.8]); plt.imshow(image, cmap='gray')
# plt.savefig('../speciale_tex/fig/clear.jpg') 
plt.figure(figsize=[4.8, 4.8]); plt.imshow(image_bin, cmap='gray')
# plt.savefig('../speciale_tex/fig/binned.jpg')



''' 
from PIL import Image
import simfun
image = np.array(Image.open('coffee.jpg').convert('L'))
image_bin = simfun.bin_sum(image, 20); plt.figure(figsize=[4.8, 4.8]); plt.imshow(image, cmap='gray')
plt.savefig('../speciale_tex/fig/clear.jpg') 
plt.figure(figsize=[4.8, 4.8]); plt.imshow(image_bin, cmap='gray'); plt.savefig('../speciale_tex/fig/binned.jpg')

for i
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


