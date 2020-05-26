#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 10:15:01 2020

@author: lasse
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

import simfun
import sys
import h5py

#CCDload = np.load('CCD3D.npy')


#a = BMB_interp(QE)
#b= BMB_interp(optics)
#c=(a+b)/2 #take mean value index-wise
#
#plt.plot(a[:,0],a[:,1],'r-')
#plt.plot(b[:,0],b[:,1],'b-')
#plt.plot(c[:,0],c[:,1],'k-')
#plt.plot(QE[0], QE[1], 'ro')
#plt.grid(1)
#plt.ylabel('Throughput')
#plt.xlabel('Wavelength $[\lambda]$')
#plt.legend(['QE', 'Optics', 'Mean', 'QE original datapoints'], loc='best')
#plt.show()
    

#QE=np.array([[350, 420, 560, 649, 789, 821, 930, 1001, 1180, 1290],[0.1, 0.7, 0.85, 0.9, 0.79, 0.85, 0.79, 0.7, 0.6, 0.3]]) 

QE = np.loadtxt('QE.txt')
#optics = np.array([[197, 200, 401, 535, 756, 831, 902, 1221],[0.99, 0.97, 0.98, 0.9, 0.99, 0.92, 0.91, 0.8]]) 
optics = np.loadtxt('optics.txt')
spec = np.loadtxt('spec2.txt', skiprows=3)

CCD_size=10240 #subpix
wl_start = 350
wl_stop = 1100
delta_lambda = 1
low_lim=300
up_lim=1150


spectrum = simfun.BMB_interp(x=spec[:,0], y=spec[:,1], wl_start=wl_start, wl_stop=wl_stop, delta_lambda=delta_lambda, kind='cubic', lowlim=low_lim, uplim=up_lim)

QE_int = simfun.BMB_interp(QE[:,0], QE[:,1], wl_start=wl_start, wl_stop=wl_stop, delta_lambda=delta_lambda, kind='cubic', lowlim=low_lim, uplim=up_lim)
optics_int = simfun.BMB_interp(optics[:,0], optics[:,1], wl_start=wl_start, wl_stop=wl_stop, delta_lambda=delta_lambda, kind='quadratic', lowlim=low_lim, uplim=up_lim)
TEC = QE_int[:,1]*optics_int[:,1]

del QE, optics, spec
#test test test

#### Generate new PSF: #########
#simfun.psf_maker(res=101, wl_start=350, wl_stop=1100, delta_lambda=1) #can also specify path and filename
################################

psf_file = h5py.File('psf_array.hdf5', 'a') #load in psf_file
psf = psf_file['psf'] #fetches the psf

'''
#Test of image generation with psf and some star position
####
#image array is WAYYY to big ~320 GB. Don't do this please
####

#file = h5py.File("image_array.hdf5", "a") #creates the file in which the arrays are stored
#file.create_dataset('image', (CCD_size, CCD_size, int((wl_stop-wl_start)/delta_lambda) ), dtype='f') #dataset for a CCD image, with stars etc.
#image = file['image']

#ech = int(np.floor(psf.shape[0]/2))
#ech2 = int(np.floor(psf.shape[1]/2))
#for i in range(100):
#    x_pos = np.random.randint(100, 10240-100)
#    y_pos = np.random.randint(100, 10240-100) #generates star position    
#    image[x_pos-ech:x_pos+ech, y_pos-ech2:y_pos+ech2, :] = image[x_pos-ech:x_pos+ech, y_pos-ech2:y_pos+ech2, :] + psf[0:100,0:100,:]
#    sys.stdout.write('.'); sys.stdout.flush(); #"Progress bar", just for visuals
#
#plt.imshow(image[:,:,-1])
'''

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




#steps=100
#ech = np.zeros((steps,2))
#for i in range(steps):
#        ech[i,0], ech[i,1] = simfun.jitter(gain=1, gain2=1)
#    
#plt.plot(ech[:,0], ech[:,1])


billede=np.zeros((1000,1000,750))
x_pos = 499
y_pos = 499
foo = psf
ux = int(np.floor(psf.shape[0]/2))
ox = int(np.floor(psf.shape[0]/2)+1)
uy = int(np.floor(psf.shape[1]/2))
oy = int(np.floor(psf.shape[1]/2)+1)
jitter = np.zeros((30,2))    
for i in range(10):
    for j in range(750):
        billede[x_pos-ux:x_pos+ox, y_pos-uy:y_pos+oy, j] = billede[x_pos-ux:x_pos+ox, y_pos-uy:y_pos+oy, j]+foo[:,:,j]
    x_j, y_j = simfun.jitter(gain=1.5 , gain2=2)
    jitter[i,0] = x_j
    jitter[i,1] = y_j
#    plt.figure()
#    plt.plot(x_j, y_j, 'r.')
    x_pos = x_pos-x_j
    y_pos = y_pos-y_j
    x_pos = int(np.around(x_pos))
    y_pos = int(np.around(y_pos))
    sys.stdout.write('*'); sys.stdout.flush(); #"Progress bar", just for visuals 
plt.figure()
plt.imshow(billede[:,:,0])

plt.figure()
plt.plot(jitter[:,0], jitter[:,1])
#### Exposure time #####










