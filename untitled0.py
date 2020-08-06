#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 10:13:21 2020

@author: lasse
"""
import numpy as np
# import simfun
# from scipy.interpolate import interp1d
# import sys 

CCDload = np.load('CCD3D.npy')
wls = np.linspace(350, 1100, 4)
# CCDload = simfun.CCD_maker(CCD_size=(100,100))
img = np.load('img.npy')
img_wl = np.load('img_wl.npy')
new_img = np.zeros((img.shape[0], img.shape[1]))
# i=0; j=0

# interp = interp1d(wls, CCDload[i,j,:])
# newCCD[i,j] = interp(img_wl[i,j])

def ccd_interp(inCCD, wls, img, img_wl):
    import sys
    from scipy.interpolate import interp1d
    
    if not wls.shape[0] is inCCD.shape[2]:
        raise TypeError("Wavelength array and input CCD depth not same size")
    if not inCCD.shape[0:2] == img.shape[0:2] == img_wl.shape:
        raise TypeError("CCD and image not same size")    
    
    for i in range(0, inCCD.shape[0]):
        for j in range(0, inCCD.shape[1]):
            interp = interp1d(wls, inCCD[i,j,:], kind="slinear", fill_value="extrapolate")
            new_img[i,j] = img[i,j]*interp(img_wl[i,j])
        sys.stdout.write('.'); sys.stdout.flush();
    return new_img

njewimg = ccd_interp(inCCD=CCDload, wls=wls, img=img, img_wl=img_wl)

'''
file_name = "tjest"
wl_endpoints=(350, 1100)
res=100
size=(100, 100)
step=1

import os
import numpy as np
# import h5py

path = os.getcwd() #Get current working directory
file_path = path + "/" + file_name +".hdf5" #Set up path to save file later

numColors = int( (wl_endpoints[1]-wl_endpoints[0])/step) # Number of colors
x_size = size[0]
y_size = size[1] #Extracts from the size input

z = np.zeros((res, res, numColors)) # Setup empty array for PSF-slices

x = np.linspace(-x_size, x_size, res) #Preparation for meshgrid
y = np.linspace(-y_size, y_size, res)
xx, yy = np.meshgrid(x, y) #define meshgrid

for i in range(wl_endpoints[0], wl_endpoints[1], step): # for-loop to create one psf for each color
    sigma_x = np.log(i)+0.5*i/100 # Used in the 2D Gaussian
    sigma_y = np.log(i)+0.5*i/100

    # 2D Gaussian function, that takes sigma_x and _y as input variables. Also takes in the meshgrid xx and yy
    zz = (1/(2*np.pi*sigma_x*sigma_y) * np.exp(-((xx)**2/(2*sigma_x**2) 
         + (yy)**2/(2*sigma_y**2))))    
    
    zz = zz/np.sum(zz) # Normalizes, so the total value (the sum of the array) =1
    z[:,:,i-350] = zz # put psf-"slice" into larger 3D array



import numpy as np
import simfun
import matplotlib.pyplot as plt
import h5py

x_size=100
y_size=100
res=100
file_name = "test_psf"
file_path = "/home/lasse/Documents/uni/Thesis/Sim/" + file_name + ".hdf5"


z = np.zeros((res, res, 750))

for i in range(350, 1100, 1):

    sigma_x = np.log(i)+0.5*i/100
    sigma_y = np.log(i)+0.5*i/100
    
    
    x = np.linspace(-x_size, x_size, res)
    y = np.linspace(-y_size, y_size, res)
    xx, yy = np.meshgrid(x, y) #define meshgrid
    
    # 2D Gaussian function, that takes sigma_x and _y as input variables. Also takes in the meshgrid xx and yy
    zz = (1/(2*np.pi*sigma_x*sigma_y) * np.exp(-((xx)**2/(2*sigma_x**2) 
         + (yy)**2/(2*sigma_y**2))))
    
    zz = zz/np.sum(zz) # Normalizes, so the total value (the sum of the array) =1
    z[:,:,i-350] = zz # put psf-"slice" into larger 3D array

# Saving the psf as a h5py file in order to store the large file
psf_file = h5py.File(file_path, "a")
psf_file.create_dataset('psf', data=z, dtype='f')

np.save("test_psf_raw.npy", z)
# plt.imshow(z)

# plt.imshow(simfun.PSF(1,1, x_size=100, y_size=100)[2])

'''


