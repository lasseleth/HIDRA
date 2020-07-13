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
    

ux = int(np.floor(psf.shape[0]/2))
ox = int(np.floor(psf.shape[0]/2)+1)
uy = int(np.floor(psf.shape[1]/2))
oy = int(np.floor(psf.shape[1]/2)+1)
#x_j, y_j = simfun.jitter(steps=1000, gain=0.9, amplitude_act=0.9, amplitude_sens=0.9)
#x_j, y_j = simfun.jitter(steps=1000, gain=0.01, amplitude_act=1.5, amplitude_sens=1.5)
#x_j, y_j = simfun.jitter(steps=10000, gain=1.9, amplitude_act=0.9, amplitude_sens=0.9)
x_j, y_j = simfun.jitter(steps=100, gain=0.2, amplitude_act=3, amplitude_sens=3)


exposure = 100 #Exposure time

''' Test of random stars in image'''
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
slitpos = [499, 499]
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

#Dispersion curve
#Straight line, steadily increasing
dispersion = -1.15*np.arange(wl_start, wl_stop, delta_lambda)+800
# plt.figure()
# plt.plot(dispersion)

eff = TEC*spectrum[:,1]

billede=np.zeros((CCD_size,CCD_size))

magni = simfun.mag(np.random.uniform(-1, 6.5)) #Random stellar brightness


# im_disp = simfun.disperser(billede, psf, magni, mask, eff, stepsize=125, plot='n') #A function that does the same as the following 
# ~60 lines
# plt.imshow(bill_disp, cmap='gray')

'''
bill_disp = np.zeros((1000,1000))
print(' ')
print('Number of colours complete:')
print('         10        20        30        40        50        60        70        80        90')

plot = str('n')

if plot == 'y':
    from matplotlib.colors import LinearSegmentedColormap
    N = 256 #8-bit value, to fix colours
    colspec = plt.cm.get_cmap('Spectral') #Fetches colourmap to use later
    vals = np.ones((N,4))
    fig1 = plt.figure()

extent= 0, 1000, 0, 1000
for k in range(0,750, 15):
    x_pos = 499 #generates star position
    y_pos = 499 #has to be here, otherwise jitter won't be applied properly
    x_0 = x_pos
    y_0 = y_pos
    for i in range(exposure):
        billede[x_pos-ux:x_pos+ox, y_pos-uy:y_pos+oy] = billede[x_pos-ux:x_pos+ox, y_pos-uy:y_pos+oy]+psf[:,:,k]*magni #adds psf values to selected area of image array
        x_pos = x_0+x_j[i]
        y_pos = y_0+y_j[i] #updates coordinates based on jitter
        x_pos = int(np.around(x_pos))
        y_pos = int(np.around(y_pos)) # rounds off the coordinates, as matrix can only take int as index
    billede_masked = billede[:,:]*mask #Overlay slit mask
    roll = np.roll(billede_masked, int(dispersion[k]), axis=1)*eff[k]
    # bill_disp = bill_disp+(np.roll(billede_masked, int(dispersion[k]), axis=1)*eff[k])+np.random.standard_normal((1000, 1000))*0.001 #Disperses the colours, using np.roll
    bill_disp = bill_disp + roll + np.random.standard_normal((1000, 1000))*0.001 #Disperses the colours, using np.roll
    
    # bill_disp = (np.roll(billede_masked, int(dispersion[k]), axis=1)*eff[k])+np.random.standard_normal((1000, 1000))*0.001 #Disperses the colours, using np.roll
    
    sys.stdout.write('/'); sys.stdout.flush(); #"Progress bar", just for visuals    
    
    ##### Plotting #####
    if plot == 'y':
        vals[:, 0] = np.linspace(0, colspec(1-k/750)[0], N) #Making new colourmap values
        vals[:, 1] = np.linspace(0, colspec(1-k/750)[1], N) #the /750 is to normalize the colormap, so values fall between 0 and 1
        vals[:, 2] = np.linspace(0, colspec(1-k/750)[2], N)
        vals[:, 3] = np.linspace(0, 1, N) #alpha, for making the cmap transparent
        newcmp = LinearSegmentedColormap.from_list(name='Spectral', colors=vals) #Creates new cmp, based on vals    
        im = plt.imshow(roll, cmap=newcmp, extent=extent) # Show array
    
    
    # vals[:, 0] = np.linspace(0, colspec(1-k/750)[0], N) #Making new colourmap values
    # vals[:, 1] = np.linspace(0, colspec(1-k/750)[1], N) #the /750 is to normalize the colormap, so values fall between 0 and 1
    # vals[:, 2] = np.linspace(0, colspec(1-k/750)[2], N)
    # vals[:, 3] = np.linspace(0, 1, N) #alpha, for making the cmap transparent
    # newcmp = LinearSegmentedColormap.from_list(name='Spectral', colors=vals) #Creates new cmp, based on vals    
    # im = plt.imshow(roll, cmap=newcmp, extent=extent) # Show array
if plot == 'y':
    plt.title('Color dispersion of sample spectrum', size=18)
    plt.xlabel('Sub-pixel', size=13)
    plt.ylabel('Sub-pixel', size=13)
    plt.savefig('disp.png', dpi=400)


'''
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

'''
jitter = simfun.jitter_im(x= x_j, y= y_j, psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]) )#Creating jitter "image"
# folded = simfun.folding(psf[:,:,0], jitter)
# plt.figure()
# # plt.plot(y_j.astype(int)+50, x_j.astype(int)+50, 'r-') #+50 to move the jitter to the correct position
# plt.imshow(folded)
# plt.imshow(jitter, alpha=0.5)
# plt.plot(y_j.mean()+50, x_j.mean()+50)
x_pos=slitpos[0]
y_pos=slitpos[1]
im_disp = np.zeros((1000,1000))
im = np.zeros((1000,1000))
from matplotlib.colors import LinearSegmentedColormap
N = 256 #8-bit value, to fix colours
colspec = plt.cm.get_cmap('Spectral') #Fetches colourmap to use later
vals = np.ones((N,4))
plot='n'
if plot=='y':
    fig1 = plt.figure()
for i in range(0,750):
    fold = simfun.folding(psf[:,:,i], jitter)
    foo = int(psf.shape[0]/2)
    im[0+x_pos-foo:len(jitter)+x_pos-foo, 0+y_pos-foo:len(jitter)+y_pos-foo] = im[0+x_pos-foo:len(jitter)+x_pos-foo, 0+y_pos-foo:len(jitter)+y_pos-foo] + fold*magni
    immask = im*mask
    roll = np.roll(immask, int(dispersion[i]), axis=1)*eff[i]
    im_disp = im_disp + roll #+ np.random.standard_normal((1000, 1000))*0.001 #Disperses the colours, using np.roll
    
    sys.stdout.write('/'); sys.stdout.flush(); #"Progress bar", just for visuals    
    
    ##### Plotting #####
    if plot == 'y':
        vals[:, 0] = np.linspace(0, colspec(i/750)[0], N) #Making new colourmap values
        vals[:, 1] = np.linspace(0, colspec(i/750)[1], N) #the /750 is to normalize the colormap, so values fall between 0 and 1
        vals[:, 2] = np.linspace(0, colspec(i/750)[2], N)
        vals[:, 3] = np.linspace(0, 1, N) #alpha, for making the cmap transparent
        newcmp = LinearSegmentedColormap.from_list(name='Spectral', colors=vals) #Creates new cmp, based on vals    
        plt.imshow(roll, cmap=newcmp) # Show array    

counts=np.array(())
for i in range(750):
    counts = np.append(counts, np.sum(im_disp[:,i]))
plt.figure(); plt.plot(counts/counts.max())
'''





mask = simfun.slit(slit_size = [slitwidth, slitheight], pos=slitpos, img_size=[10240, 10240])

jitter = simfun.jitter_im(x= x_j, y= y_j, psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]) )#Creating jitter "image"

def disperser2(jitter, psf, pos, image_size, dispersion, eff, magni, steps=50, plot='n'):
    # jitter = simfun.jitter_im(x= x_j, y= y_j, psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]) )#Creating jitter "image"
    x_pos=pos[0]
    y_pos=pos[1]
    im_disp = np.zeros((10240,10240))
    im = np.zeros((10240,10240))
    from matplotlib.colors import LinearSegmentedColormap
    N = 256 #8-bit value, to fix colours
    colspec = plt.cm.get_cmap('Spectral') #Fetches colourmap to use later
    vals = np.ones((N,4))
    # plot='n'
    if plot=='y':
        plt.figure()
    for i in range(0,750, steps):
        fold = simfun.folding(psf[:,:,i], jitter)
        foo = int(psf.shape[0]/2)
        im[0+x_pos-foo:len(jitter)+x_pos-foo, 0+y_pos-foo:len(jitter)+y_pos-foo] = im[0+x_pos-foo:len(jitter)+x_pos-foo, 0+y_pos-foo:len(jitter)+y_pos-foo] + fold*magni
        immask = im*mask
        roll = np.roll(immask, int(dispersion[i]), axis=1)*eff[i]
        im_disp = im_disp + roll #+ np.random.standard_normal((1000, 1000))*0.001 #Disperses the colours, using np.roll
        
        sys.stdout.write('/'); sys.stdout.flush(); #"Progress bar", just for visuals    
        
        ##### Plotting #####
        if plot == 'y':
            vals[:, 0] = np.linspace(0, colspec(i/750)[0], N) #Making new colourmap values
            vals[:, 1] = np.linspace(0, colspec(i/750)[1], N) #the /750 is to normalize the colormap, so values fall between 0 and 1
            vals[:, 2] = np.linspace(0, colspec(i/750)[2], N)
            vals[:, 3] = np.linspace(0, 1, N) #alpha, for making the cmap transparent
            newcmp = LinearSegmentedColormap.from_list(name='Spectral', colors=vals) #Creates new cmp, based on vals    
            plt.imshow(roll, cmap=newcmp) # Show array    
    
    if plot=='y':
        plt.title('Color dispersion of sample spectrum', size=18)
        plt.xlabel('Sub-pixel', size=13)
        plt.ylabel('Sub-pixel', size=13)
    counts=np.array(())
    for i in range(750):
        counts = np.append(counts, np.sum(im_disp[:,i]))
    if plot=='y': plt.figure(); plt.plot(counts/counts.max())

    return im_disp

# test=disperser2(jitter, psf, pos=slitpos, image_size=(1000,1000), dispersion=dispersion, eff=eff, magni=magni, steps=125, plot='y')
















