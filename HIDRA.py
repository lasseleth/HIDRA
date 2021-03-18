'''
Hyperspectral Instrument Data Resemblance Algorithm - A simple modular Python simulator for low-resolution spectroscopy
For information about inputs, installation etc., please see thesis work, chapter 7.
Made by Lasse L. S. Berthelsen, as a part of the thesis work. 
2020

'''


import numpy as np
import matplotlib.pyplot as plt
import input_file as inp
import simfun
import scipy.signal

def rms(x):
    return np.sqrt(np.mean(x**2))


def int_r(r1, r2, rang):
    from scipy.interpolate import interp1d
    x= np.arange(len(r1))
    xnew = np.arange(0, len(r1), 0.001)
    f1 = interp1d(x, r1, kind=3, fill_value="extrapolate")
    f2 = interp1d(x, r2, kind=3, fill_value="extrapolate")
    
    r1_int = f1(xnew)
    r2_int = f2(xnew)
    return r1_int, r2_int
#### SETUP PHASE ####

# del x_j, y_j


spec_eff, spec_eff2, jitter, x_j, y_j, psf, img_size, sub_pixel, pl_arc_mm, disper, mask, slitpos, background = simfun.setup(inp)

in_spec = np.loadtxt(inp.in_spec)
in_spec2 = np.loadtxt(inp.in_spec2)
wl_ran = inp.wl_ran
exp = inp.exp
slit = inp.slit
CCD = np.load(inp.in_CCD)

image, image_wl=simfun.disperser(wl_endpoints=wl_ran, jit_img=jitter, psf_ends=[15, 45], pos=slitpos, image_size=img_size, 
                                        dispersion=disper, eff=spec_eff, mask_img=mask, steps=1, plot='n')
# New jitter files:
# for i in range(1, 21):
    # x_j, y_j = simfun.func_jitter(entries=(exp*inp.step), gain=0.15, dt=5)
    # j = np.vstack([x_j, y_j])
    # np.save("runs/jitter/j"+str(i), j)


#### SETUP PHASE COMPLETE ####
#### IMAGE FORMATION BEGINS ####
def the_thing(image, image2, sub_pixel=sub_pixel, wl_ran=wl_ran, disper=disper, slitpos=slitpos, move="y", noiseinp="n"):
    rout = simfun.bin_sum(image, sub_pixel)
    r1 = simfun.read_out(rout) 
    
    rin = simfun.bin_sum(image2, sub_pixel)
    r2 = simfun.read_out(rin) #sum image and readout 
    r1, r2 = int_r(r1, r2, wl_ran) #Interpolate to higher resolution
    print("\nBinning and read-out done")
    
    if noiseinp == "y":
        #rs= np.random.RandomState()
        # no = rs.poisson(np.mean(image)*sub_pixel*sub_pixel,  size=(int(image.shape[0] /sub_pixel), int(image.shape[1]/sub_pixel)))
        # ni = rs.poisson(np.mean(image2)*sub_pixel*sub_pixel, size=(int(image2.shape[0]/sub_pixel), int(image2.shape[1]/sub_pixel)))
        no =  noise1d(rout)
        ni = noise1d(rin)
        
        
        
    if move == "y":
        autocor = scipy.signal.correlate(r1, r1, mode="same") #perform autocorr.
        cor = scipy.signal.correlate(r1, r2, mode="same") #Regular correlation
        first = np.argmax(autocor)
        second = np.argmax(cor)
        delta = first-second #amount of sub-pixels to move r1, for the two spectra to overlap
        
        if noiseinp == "y":
            rout = simfun.read_out(simfun.bin_sum(image*CCD, sub_pixel)+no) 
            rin = simfun.read_out(simfun.bin_sum(image2*CCD, sub_pixel)+ni)
            r1, r2 = int_r(rout, rin, wl_ran)
        
        r1 = np.roll(r1, delta) #Move r1
        
        del first, second, autocor, cor, rout, rin
    if not move == "y":
        delta = 0
    # noiseinp = input("Include photon noise? y/n: ")
    # noiseinp = "n"    
    
    # r1 = scipy.ndimage.filters.uniform_filter1d(r1, size=5000)
    # r2 = scipy.ndimage.filters.uniform_filter1d(r2, size=5000)
    
    pos = (disper[0]+slitpos[0])*100.0 #Position of each wavelength on the detector
    from scipy.stats import linregress
    a, b, r, p, s = linregress(pos, np.arange(wl_ran[0], wl_ran[1])) #Linear regression to find the lambda/pixel correlation
    wavelength = a*np.arange(img_size[1]*100.0)+(b)
    del a, b, r, p, s, 
    
    wave = wavelength[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1] #remove outlying entries, where the spectrum is not present (belo 300 nm, and above 1000) 
    r1 = r1[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1] 
    r2 = r2[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1]
    # plt.plot(in_spec[:,0], in_spec[:,1]/in_spec2[:,1])
    print("\nMoving mean filters: ")
    from astropy.convolution import convolve, Gaussian1DKernel
    r1 = convolve(r1,kernel = Gaussian1DKernel(4246.6)) #Moving Mean filter by convolution. Kernel is Gaussian, input is sigma
    # print("Done! \nFirst spectrum complete, continuing to second:\n")
    r2 = convolve(r2,kernel = Gaussian1DKernel(4246.6)) #https://docs.astropy.org/en/stable/convolution/
    print("MMF done\n")
    plt.plot(wave, (r1-r2)/r1)
    return r1, r2, wave, delta



# ======================================================================================================================
# slit = ['pix', 2, 3.5]
# slit_size = simfun.convert_slit(unit = slit[0], size = slit[1:3], convert_factor = pl_arc_mm) #Convert slit size to pixels
# slit_size[0] = slit_size[0]*sub_pixel #Convert to subpixels
# slit_size[1] = slit_size[1]*sub_pixel
# slitpos = [150, 249] #Slit position on the sub-pixel CCD image. Arbitrary position..
# mask = simfun.func_slit(slit_size = np.floor(slit_size).astype(int), pos=slitpos, image_size=img_size)

# mask = np.ones((img_size))
# for i in range(1, 11):
#     jitter = np.load("jitter/j"+str(i)+".npy")
#     jitter = simfun.jitter_im(x= jitter[0,:], y= jitter[1,:], psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]) )
    
#     jitter2 = np.load("jitter/j"+str(i+10)+".npy")
#     jitter2 = simfun.jitter_im(x= jitter2[0,:], y= jitter2[1,:], psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]) )
    
    
    
#     image1, image1_wl=simfun.disperser(wl_endpoints=wl_ran, jit_img=jitter, psf_img=psf, pos=slitpos, image_size=img_size, 
#                                       dispersion=disper, eff=spec_eff, mask_img=mask, steps=1, plot='n')
    
#     image2, image2_wl=simfun.disperser(wl_endpoints=wl_ran, jit_img=jitter2, psf_img=psf, pos=slitpos, image_size=img_size, 
#                                         dispersion=disper, eff=spec_eff2, mask_img=mask, steps=1, plot='n')
    
#     #np.save("e600/slit"+str(slit[1])+"/r"+str(i)+"/out", image1)
#     np.save("e300_bigpsf_noslit/r"+str(i)+"/out", image1)
#     np.save("e300_bigpsf_noslit/r"+str(i)+"/in", image2)

# # ======================================================================================================================
# slit = ['pix', 2, 3.5]
# slit_size = simfun.convert_slit(unit = slit[0], size = slit[1:3], convert_factor = pl_arc_mm) #Convert slit size to pixels
# slit_size[0] = slit_size[0]*sub_pixel #Convert to subpixels
# slit_size[1] = slit_size[1]*sub_pixel
# slitpos = [150, 249] #Slit position on the sub-pixel CCD image. Arbitrary position..
# mask = simfun.func_slit(slit_size = np.floor(slit_size).astype(int), pos=slitpos, image_size=img_size)


# for i in range(1, 11):
#     jitter = np.load("jitter/j"+str(i)+".npy")
#     jitter = simfun.jitter_im(x= jitter[:,0], y= jitter[:,1], psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]) )
    
#     jitter2 = np.load("jitter/j"+str(i+10)+".npy")
#     jitter2 = simfun.jitter_im(x= jitter2[:,0], y= jitter2[:,1], psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]) )
    
    
#     image1, image1_wl=simfun.disperser(wl_endpoints=wl_ran, jit_img=jitter, psf_img=psf, pos=slitpos, image_size=img_size, 
#                                       dispersion=disper, eff=spec_eff, mask_img=mask, steps=1, plot='n')
    
#     image2, image2_wl=simfun.disperser(wl_endpoints=wl_ran, jit_img=jitter2, psf_img=psf, pos=slitpos, image_size=img_size, 
#                                         dispersion=disper, eff=spec_eff2, mask_img=mask, steps=1, plot='n')
    
#     np.save("1200/slit"+str(slit[1])+"/r"+str(i)+"/out", image1)
#     np.save("1200/slit"+str(slit[1])+"/r"+str(i)+"/in", image2)

# # ======================================================================================================================
# slit = ['pix', 4, 3.5]
# slit_size = simfun.convert_slit(unit = slit[0], size = slit[1:3], convert_factor = pl_arc_mm) #Convert slit size to pixels
# slit_size[0] = slit_size[0]*sub_pixel #Convert to subpixels
# slit_size[1] = slit_size[1]*sub_pixel
# slitpos = [150, 249] #Slit position on the sub-pixel CCD image. Arbitrary position..
# mask = simfun.func_slit(slit_size = np.floor(slit_size).astype(int), pos=slitpos, image_size=img_size)

# for i in range(1, 11):
#     jitter = np.load("runs/jitter/j"+str(i)+".npy")
#     jitter = simfun.jitter_im(x= jitter[:,0], y= jitter[:,1], psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]) )
    
#     jitter2 = np.load("runs/jitter/j"+str(i+10)+".npy")
#     jitter2 = simfun.jitter_im(x= jitter2[:,0], y= jitter2[:,1], psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]) )
    
    
#     image1, image1_wl=simfun.disperser(wl_endpoints=wl_ran, jit_img=jitter, psf_img=psf, pos=slitpos, image_size=img_size, 
#                                      dispersion=disper, eff=spec_eff, mask_img=mask, steps=1, plot='n')
    
#     image2, image2_wl=simfun.disperser(wl_endpoints=wl_ran, jit_img=jitter2, psf_img=psf, pos=slitpos, image_size=img_size, 
#                                        dispersion=disper, eff=spec_eff2, mask_img=mask, steps=1, plot='n')
    
#     np.save("runs/slit"+str(slit[1])+"/r"+str(i)+"/out", image1)
#     np.save("runs/slit"+str(slit[1])+"/r"+str(i)+"/in", image2)

# # ======================================================================================================================
# slit = ['pix', 12, 12]
# slit_size = simfun.convert_slit(unit = slit[0], size = slit[1:3], convert_factor = pl_arc_mm) #Convert slit size to pixels
# slit_size[0] = slit_size[0]*sub_pixel #Convert to subpixels
# slit_size[1] = slit_size[1]*sub_pixel
# slitpos = [150, 249] #Slit position on the sub-pixel CCD image. Arbitrary position..
# mask = simfun.func_slit(slit_size = np.floor(slit_size).astype(int), pos=slitpos, image_size=img_size)

# for i in range(1, 11):
#     jitter = np.load("runs/jitter/j"+str(i)+".npy")
#     jitter = simfun.jitter_im(x= jitter[:,0], y= jitter[:,1], psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]) )
    
#     jitter2 = np.load("runs/jitter/j"+str(i+10)+".npy")
#     jitter2 = simfun.jitter_im(x= jitter2[:,0], y= jitter2[:,1], psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]) )
    
    
#     image1, image1_wl=simfun.disperser(wl_endpoints=wl_ran, jit_img=jitter, psf_img=psf, pos=slitpos, image_size=img_size, 
#                                      dispersion=disper, eff=spec_eff, mask_img=mask, steps=1, plot='n')
    
#     image2, image2_wl=simfun.disperser(wl_endpoints=wl_ran, jit_img=jitter2, psf_img=psf, pos=slitpos, image_size=img_size, 
#                                        dispersion=disper, eff=spec_eff2, mask_img=mask, steps=1, plot='n')
    
#     np.save("runs/slit"+str(slit[1])+"/r"+str(i)+"/out", image1)
#     np.save("runs/slit"+str(slit[1])+"/r"+str(i)+"/in", image2)


'''
check = input("Run simulation? y/n: ")

if check == "y":
    image, image_wl=simfun.disperser(wl_endpoints=wl_ran, jit_img=jitter, psf_img=psf, pos=slitpos, image_size=img_size, 
                                        dispersion=disper, eff=spec_eff, mask_img=mask, steps=1, plot='n')

    # lam = (image_wl+np.mean([wl_ran[0], wl_ran[1]])*10**(-80))/(image+10**(-80)) #To find the dispersion wavelength (for calibrating later)
    # summed = simfun.bin_sum(image, sub_pixel)
    # r1 = simfun.read_out(summed)

    number = round(rms(y_j*(pl_arc_mm/sub_pixel)), 3)
    if "." in str(number):
        jin = str(number).split(".")[0]
        jid = str(number).split(".")[-1]
    
    number = round(slit[1]*(pl_arc_mm), 3)
    if "." in str(number):
        sn = str(number).split(".")[0]
        sd = str(number).split(".")[-1]
    name = 'out'+'-e'+str(exp)+'-j'+jin+'_'+jid+'-s'+sn+'_'+sd
    
    np.save('runs/'+name, image)1
    np.save('runs/'+name+'_wl', image_wl)


    #New jitter, and redo simulation, now with spec_eff2
    x_j2, y_j2 = simfun.func_jitter(entries=(exp*inp.step), gain=0.15, dt=5)
    # x_j2, y_j2 = simfun.jitter(entries=(exp*step), gain=0.02, dt=10)
    jitter2 = np.stack((x_j2, y_j2), axis=-1)
    jitter2 = simfun.jitter_im(x= jitter2[:,0], y= jitter2[:,1], psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]) )
    
    image2, image2_wl=simfun.disperser(wl_endpoints=wl_ran, jit_img=jitter2, psf_img=psf, pos=slitpos, image_size=img_size, 
                                            dispersion=disper, eff=spec_eff2, mask_img=mask, steps=1, plot='n')
    
    number = round(rms(y_j2*(pl_arc_mm/sub_pixel)), 3)
    if "." in str(number):
        jin = str(number).split(".")[0]
        jid = str(number).split(".")[-1]
    
    number = round(slit[1]*(pl_arc_mm), 3)
    if "." in str(number):
        sn = str(number).split(".")[0]
        sd = str(number).split(".")[-1]
    name2 = 'in'+'-e'+str(exp)+'-j'+jin+'_'+jid+'-s'+sn+'_'+sd
    
    np.save('runs/'+name2, image2)
    np.save('runs/'+name2+'_wl', image2_wl)

'''


# control_in = np.load('runs/in-e300-j1_357-s7_054.npy')
# sumcin = simfun.bin_sum(control_in, 10)
# rcin = simfun.read_out(sumcin)

# control_out = np.load('runs/out-e300-j1_372-s7_054.npy')
# sumcout = simfun.bin_sum(control_out, 10)
# rcout = simfun.read_out(sumcout)

# jit_out = np.load('runs/out-e300-j2_407-s7_054.npy')
# sumjitout = simfun.bin_sum(jit_out, 10)
# rjitout = simfun.read_out(sumjitout)

# jit_in = np.load('runs/in-e300-j2_853-s7_054.npy')
# sumjitin = simfun.bin_sum(jit_in, 10)
# rjitin = simfun.read_out(sumjitin)

# rcout = rjitout
# rcin = rjitin

# from scipy.signal import find_peaks
# peak_ref, _ = find_peaks(in_spec[:,1]*(-1), threshold=4, prominence=200, distance=3)
# # peaks, _ = find_peaks(-1*rcout, threshold=2.6e5, distance=3, prominence=500000)
# peaks, _ = find_peaks(-1*rcout, threshold=1.65e5, distance=1, prominence=500000)
# peaks = peaks[:-1]
# peak_ref = peak_ref[1:]
# from scipy.stats import linregress #Function to find a and b for a linear dispersion function. This of course only works as the dispersion is linear!
# a, b, r, p, s = linregress(peaks, in_spec[:,0][peak_ref])
# wavelength = a*np.arange(img_size[1]/10)+(b) #Making the new wavelength array for the readout values.
# del a, b, r, p, s
# wave = wavelength[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1]

# r1 = rcout[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1]
# r2 = rcin[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1]

# bin_width=100
# cr1 = [sum(r1[i:i+bin_width]) for i in range(50, len(r1), bin_width)]
# cr2 = [sum(r2[i:i+bin_width]) for i in range(50, len(r2), bin_width)]

# foo = np.zeros((int(r1.shape[0]/bin_width), 9))
# for i in range(7):
#     a = cr1[i]
#     foo[i] = np.array([cr1[i]/a, cr1[i]/a, cr1[i]/a, cr2[i]/a, cr2[i]/a, cr2[i]/a, cr1[i]/a, cr1[i]/a, cr1[i]/a ])

# color = ['purple', 'blue', 'cyan', 'green', 'yellow', 'orange', 'red']
# for i in range(7): plt.plot(foo[i,:], '.', markersize=13, color=color[i])
# for i in range(7): plt.plot(foo[i,:], '-', color=color[i])
# plt.xlabel("Time, arbitrary units", size=14)
# plt.ylabel("Relative luminosity", size=14)
# plt.tick_params(labelsize=14)
# plt.grid()
def noise1d(x, RON=5):
    noise = np.zeros((x.shape))
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            noise[i,j] = (np.sqrt(x[i,j])+RON)*np.random.normal(0,1)
    return noise

'''
def the_thing(image, image2, sub_pixel=sub_pixel, wl_ran=wl_ran, disper=disper, slitpos=slitpos, move="y", noiseinp="n"):
    rout = simfun.read_out(simfun.bin_sum(image, sub_pixel))  
    rin = simfun.read_out(simfun.bin_sum(image2, sub_pixel)) #sum image and readout 
    r1, r2 = int_r(rout, rin, wl_ran) #Interpolate to higher resolution
    print("\nBinning and read-out done")
    if move == "y":
        autocor = scipy.signal.correlate(r1, r1, mode="same") #perform autocorr.
        cor = scipy.signal.correlate(r1, r2, mode="same") #Regular correlation
        first = np.argmax(autocor)
        second = np.argmax(cor)
        delta = first-second #amount of sub-pixels to move r1, for the two spectra to overlap
        
        r1 = np.roll(r1, delta) #Move r1
        
        
        
        del first, second, autocor, cor
    if not move == "y":
        delta = 0
    # noiseinp = input("Include photon noise? y/n: ")
    # noiseinp = "n"
    if noiseinp == "y":
        rs= np.random.RandomState()
        no = rs.poisson(np.mean(image),  size=(int(image.shape[0] /sub_pixel), 700))
        ni = rs.poisson(np.mean(image2), size=(int(image2.shape[0]/sub_pixel), 700))
    
        no = simfun.read_out(no)
        ni = simfun.read_out(ni)
    
    
    # r1 = scipy.ndimage.filters.uniform_filter1d(r1, size=5000)
    # r2 = scipy.ndimage.filters.uniform_filter1d(r2, size=5000)
    
    pos = (disper[0]+slitpos[0])*100.0 #Position of each wavelength on the detector
    from scipy.stats import linregress
    a, b, r, p, s = linregress(pos, np.arange(wl_ran[0], wl_ran[1])) #Linear regression to find the lambda/pixel correlation
    wavelength = a*np.arange(img_size[1]*100.0)+(b)
    del a, b, r, p, s, 
    
    wave = wavelength[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1] #remove outlying entries, where the spectrum is not present (belo 300 nm, and above 1000) 
    r1 = r1[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1] 
    r2 = r2[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1]
    # plt.plot(in_spec[:,0], in_spec[:,1]/in_spec2[:,1])
    print("\nMoving mean filters: ")
    from astropy.convolution import convolve, Gaussian1DKernel
    r1 = convolve(r1,kernel = Gaussian1DKernel(4246.6)) #Moving Mean filter by convolution. Kernel is Gaussian, input is sigma
    # print("Done! \nFirst spectrum complete, continuing to second:\n")
    r2 = convolve(r2,kernel = Gaussian1DKernel(4246.6)) #https://docs.astropy.org/en/stable/convolution/
    print("MMF done\n")
    plt.plot(wave, r1/r2)
    return r1, r2, wave, delta
'''


# r5out, r5in, wave = the_thing(image, image2)

# ri = np.load('slit2/r'+str(1)+'/in.npy')
# ro = np.load('slit2/r'+str(1)+'/out.npy')
# rout = simfun.read_out(simfun.bin_sum(ro, sub_pixel))  
# rin = simfun.read_out(simfun.bin_sum(ri, sub_pixel))
# pos = (disper[0]+slitpos[0]) #Position of each wavelength on the detector
# from scipy.stats import linregress
# a, b, r, p, s = linregress(pos, np.arange(wl_ran[0], wl_ran[1])) #Linear regression to find the lambda/pixel correlation
# wavelength = a*np.arange(img_size[1])+(b)
# wavelength = wavelength[0::10]
# wave = wavelength[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1] #remove outlying entries, where the spectrum is not present (belo 300 nm, and above 1000) 
# r1 = rout[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1] 
# r2 = rin[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1]
# print("interpolating")
# r1, r2 = int_r(r1, r2, wl_ran) #Interpolate to higher resolution
# print(" done")
# autocor = scipy.signal.correlate(r1, r1, mode="same") #perform autocorr.
# cor = scipy.signal.correlate(r1, r2, mode="same") #Regular correlation
# first = np.argmax(autocor)
# second = np.argmax(cor)
# delta = first-second #amount of sub-pixels to move r1, for the two spectra to overlap
# print("MMF")
# r1 = np.roll(r1, delta) #Move r1
# from astropy.convolution import convolve, Gaussian1DKernel
# r1 = convolve(r1,kernel = Gaussian1DKernel(4246)) #Moving Mean filter by convolution. Kernel is Gaussian, input is sigma
# # print("Done! \nFirst spectrum complete, continuing to second:\n")
# r2 = convolve(r2,kernel = Gaussian1DKernel(4246)) #https://docs.astropy.org/en/stable/convolution/



# r1o = np.load('runs/s2/out-e300-j1_33-s9_406.npy')
# r1i = np.load('runs/s2/in-e300-j1_333-s9_406.npy')
# n1o = simfun.noise(size=r1o.shape, image=r1o)
# n1i = simfun.noise(size=r1i.shape, image=r1i)
# r1i = r1i+n1i
# r1o = r1o+n1o
# del n1o, n1i
# r1o, r1i, wave = the_thing(r1o, r1i)



# r2i = np.load('runs/s2/in-e300-j1_323-s9_406.npy')
# r2o = np.load('runs/s2/out-e300-j1_369-s9_406.npy')
# n2o = simfun.noise(size=r2o.shape, image=r2o)
# n2i = simfun.noise(size=r2i.shape, image=r2i)
# r2i = r2i+n2i
# r2o = r2o+n2o
# del n2o, n2i
# r2o, r2i, wave = the_thing(r2o, r2i)


# r3i = np.load('runs/s2/in-e300-j1_41-s9_406.npy')
# r3o = np.load('runs/s2/out-e300-j1_368-s9_406.npy')
# n3o = simfun.noise(size=r3o.shape, image=r3o)
# n3i = simfun.noise(size=r3i.shape, image=r3i)
# r3i = r3i+n3i
# r3o = r3o+n3o
# del n3o, n3i
# r3o, r3i, wave = the_thing(r3o, r3i)


# ri_s2 = np.vstack([r1i, r2i, r3i])
# ro_s2 = np.vstack([r1o, r2o, r3o])
# del r1i, r2i, r3i
# del r1o, r2o, r3o


# ri_s2_std = np.zeros((ri_s2.shape[1]))
# ro_s2_std = np.zeros((ro_s2.shape[1]))
# for i in range(len(ri_s2_std)):
#     ri_s2_std[i] = np.std(ri_s2[:,i])
#     ro_s2_std[i] = np.std(ro_s2[:,i])


# ro_s2_mean = np.zeros((ro_s2.shape[1]))
# ri_s2_mean = np.zeros((ri_s2.shape[1]))
# for i in range(ri_s2.shape[1]):
#     ri_s2_mean[i] = np.mean(ri_s2[:,i])
#     ro_s2_mean[i] = np.mean(ro_s2[:,i])

# plt.figure()
# plt.plot(wave, ro_s2_mean/ri_s2_mean)
# plt.plot(in_spec[:,0], in_spec[:,1]/in_spec2[:,1])
# # plt.fill_between(wave, ((ro_mean-ro_std)/(ri_mean-ri_std)), ((ro_mean+ro_std)/(ri_mean+ri_std)), alpha=0.2)
# plt.legend(["Output", "Input", "1 $\sigma$"], fontsize=14)
# plt.tick_params(labelsize=12)
# plt.xlabel("Wavelength [$nm$]", size=14)
# plt.ylabel("Ratio", size=14)
# plt.grid()

'''

r1i = np.load('runs/s1/in-e300-j1_346-s4_703.npy')
r1o = np.load('runs/s1/out-e300-j1_375-s4_703.npy')
n1o = simfun.noise(size=r1o.shape, image=r1o)
n1i = simfun.noise(size=r1i.shape, image=r1i)
r1i = r1i*CCD+n1i
r1o = r1o*CCD+n1o
del n1o, n1i
r1o, r1i, wave = the_thing(r1o, r1i)



r2i = np.load('runs/s1/in-e300-j1_37-s4_703.npy')
r2o = np.load('runs/s1/out-e300-j1_376-s4_703.npy')
n2o = simfun.noise(size=r2o.shape, image=r2o)
n2i = simfun.noise(size=r2i.shape, image=r2i)
r2i = r2i*CCD+n2i
r2o = r2o*CCD+n2o
del n2o, n2i
r2o, r2i, wave = the_thing(r2o, r2i)


r3i = np.load('runs/s1/in-e300-j1_394-s4_703.npy')
r3o = np.load('runs/s1/out-e300-j1_395-s4_703.npy')
n3o = simfun.noise(size=r3o.shape, image=r3o)
n3i = simfun.noise(size=r3i.shape, image=r3i)
r3i = r3i*CCD+n3i
r3o = r3o*CCD+n3o
del n3o, n3i
r3o, r3i, wave = the_thing(r3o, r3i)


r4i = np.load('runs/s1/in-e300-j1_399-s4_703.npy')
r4o = np.load('runs/s1/out-e300-j1_409-s4_703.npy')
n4o = simfun.noise(size=r4o.shape, image=r4o)
n4i = simfun.noise(size=r4i.shape, image=r4i)
r4i = r4i*CCD+n4i
r4o = r4o*CCD+n4o
del n4o, n4i
r4o, r4i, wave = the_thing(r4o, r4i)


ri = np.vstack([r1i, r2i, r3i, r4i])
ro = np.vstack([r1o, r2o, r3o, r4o])
del r1i, r2i, r3i, r4i
del r1o, r2o, r3o, r4o


ri_std = np.zeros((ri.shape[1]))
ro_std = np.zeros((ro.shape[1]))
for i in range(len(ri_std)):
    ri_std[i] = np.std(ri[:,i])
    ro_std[i] = np.std(ro[:,i])


ro_mean = np.zeros((ro.shape[1]))
ri_mean = np.zeros((ri.shape[1]))
for i in range(ri.shape[1]):
    ri_mean[i] = np.mean(ri[:,i])
    ro_mean[i] = np.mean(ro[:,i])

cm = 1/2.54
plt.figure(figsize=(19.5*cm, 13.8*cm))
plt.plot(wave, ro_mean/ri_mean)
plt.plot(in_spec[:,0], in_spec[:,1]/in_spec2[:,1])
plt.fill_between(wave, ((ro_mean-ro_std)/(ri_mean-ri_std)), ((ro_mean+ro_std)/(ri_mean+ri_std)), alpha=0.2)
plt.legend(["Output", "Input", "1 $\sigma$"], fontsize=12, loc="lower center")
plt.tick_params(labelsize=12)
plt.xlabel("Wavelength [$nm$]", size=12)
plt.ylabel("Ratio", size=12)
plt.grid()
plt.tight_layout()
#plt.savefig("../speciale_tex/fig/noise2.png", dpi=300)













=========================================
r1, r2 = int_r(rcout, rcin)

autocor = scipy.signal.correlate(r1, r1, mode="same")
cor = scipy.signal.correlate(r1, r2, mode="same")
first = np.argmax(autocor)
second = np.argmax(cor)
delta = first-second
r1 = np.roll(r1, delta)

from astropy.convolution import convolve, Gaussian1DKernel
r1 = convolve(r1,kernel = Gaussian1DKernel(5000))
r2 = convolve(r2,kernel = Gaussian1DKernel(5000)) #Moving Mean filter by convolution. 
# r1 = scipy.ndimage.filters.uniform_filter1d(r1, size=5000)
# r2 = scipy.ndimage.filters.uniform_filter1d(r2, size=5000)

pos = (disper[0]+slitpos[0])*100
from scipy.stats import linregress
a, b, r, p, s = linregress(pos, np.arange(wl_ran[0], wl_ran[1]))
wavelength = a*np.arange(img_size[1]*100)+(b)

wave = wavelength[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1]
r1 = r1[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1]
r2 = r2[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1]
plt.plot(in_spec[:,0], in_spec[:,1]/in_spec2[:,1])
plt.plot(wave, r1/r2)
=========================================


control_in = np.load('runs/in-e300-j1_357-s7_054.npy')
sumcin = simfun.bin_sum(control_in, 10)
rcin = simfun.read_out(sumcin)

control_out = np.load('runs/out-e300-j1_372-s7_054.npy')
sumcout = simfun.bin_sum(control_out, 10)
rcout = simfun.read_out(sumcout)

exp_out = np.load('runs/out-e600-j1_403-s7_054.npy')
sumexpout = simfun.bin_sum(exp_out, 10)
rexpout = simfun.read_out(sumexpout)

exp_in = np.load('runs/in-e600-j1_464-s7_054.npy')
sumexpin = simfun.bin_sum(exp_in, 10)
rexpin = simfun.read_out(sumexpin)

slit_in = np.load('runs/in-e300-j1_4-s49_381.npy')
sumslitin = simfun.bin_sum(slit_in, 10)
rslitin = simfun.read_out(sumslitin)

slit_out = np.load('runs/out-e300-j1_383-s49_381.npy')
sumslitout = simfun.bin_sum(slit_out, 10)
rslitout = simfun.read_out(sumslitout)

jit_out = np.load('runs/out-e300-j2_407-s7_054.npy')
sumjitout = simfun.bin_sum(jit_out, 10)
rjitout = simfun.read_out(sumjitout)

jit_in = np.load('runs/in-e300-j2_853-s7_054.npy')
sumjitin = simfun.bin_sum(jit_in, 10)
rjitin = simfun.read_out(sumjitin)



lam2 = (image2_wl+np.mean([wl_ran[0], wl_ran[1]])*10**(-80))/(image2+10**(-80)) #To find the dispersion wavelength (for calibrating later)

summed2 = simfun.bin_sum(image2, sub_pixel)
r2 = simfun.read_out(summed2)


from scipy.signal import find_peaks # load function to find spectral lines
peak_ref, _ = find_peaks(in_spec[:,1]*(-1), threshold=4, prominence=200, distance=3) #Finds most prominent peaks. Values found empirically (trial and error)

# rr = simfun.read_out(image) 
peaks, _ = find_peaks(-1*r1,threshold=2.6e5, distance=3, prominence=500000) #Peaks of the readout. Again, trial and error with the values

from scipy.stats import linregress #Function to find a and b for a linear dispersion function. This of course only works as the dispersion is linear!
a, b, r, p, s = linregress(peaks, in_spec[:,0][peak_ref[1:]]) 
wavelength = a*np.arange(img_size[1]/10)+(b) #Making the new wavelength array for the readout values.
del a, b, r, p, s

r1 = r1[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-2]
wave = wavelength[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-2]

r2 = r2[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-2]

'''



"""
===============================================================
                        Plot bg_spec.png
===============================================================
bg_spec = np.loadtxt(inp.bg_spec)
if not bg_spec.shape[0] == wl_ran[1]-wl_ran[0]: #interpolate if values are missing
    bg_spec = simfun.interp(x=bg_spec[:,0], y=bg_spec[:,1], wl_ran=wl_ran, kind='cubic', lowlim=wl_ran[0]-50, uplim=wl_ran[1]+50)
detector_area = (pl_arc_mm*img_size[0]/sub_pixel)*(pl_arc_mm*img_size[1]/sub_pixel) #Collecting area of the detector measured in arcsec^2
bg_spec[:,1] = bg_spec[:,1]*detector_area #Multiply by detector area
plt.figure()
plt.plot(bg_spec[:,0], bg_spec[:,1], linewidth=4)
plt.grid()
plt.xlabel("Wavelength [$nm$]", size=12)
plt.ylabel("Intensity [$photons$ $s^{-1}$ $nm^{-1}$]", size=12)
plt.tick_params(labelsize=12)
plt.tight_layout()

# plt.plot(wave, r1)
# plt.plot(wave, r2)
# plt.plot(in_spec[:,0], in_spec[:,1]*50000)
===============================================================



===============================================================
                    Plot eta_in.png
===============================================================
foo = 1
eta_in = inp.eta_in
plt.figure()
args = eta_in
for x in args: #takes several input arrays
    loaded = np.loadtxt(x)
    if not loaded.shape[0] == wl_ran[1]-wl_ran[0]: #if input is not of the correct length, this will interpolate
        temp = simfun.interp(loaded[:,0], loaded[:,1], wl_ran=wl_ran, kind='cubic', lowlim=wl_ran[0]-50, uplim=wl_ran[1]+50)[:,1]
    else:
        temp = loaded
    foo = foo * temp
    plt.plot(np.arange(300, 1000), temp, linewidth=4)
eta_in=foo
plt.plot(np.arange(300,1000), eta_in, linewidth=4)
plt.legend(["CCD QE", "$\eta$ of telescope optics", "Total Efficiency Curve"], loc="lower center", fontsize=12)
plt.grid()
plt.xlabel("Wavelength [nm]", size=12)
plt.ylabel("Relative spectral throughput, $\eta$", size=12)
===============================================================



===============================================================
                    Plot ratio_s4_vs_s12.png
===============================================================
fig = plt.figure(figsize=(9,5))
ax1  = fig.add_subplot(211)
ax1.plot(wave, (roall4-riall4)/roall4)
plt.axis([ 265.9, 1035, -0.00333, 0.02540])
plt.ylabel("Normalized ratio", size=13)
plt.tick_params(labelbottom=False, labelsize=12)
plt.grid()

color=['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'tab:brown']
ax2  = fig.add_subplot(212)
for i in range(10):
    ax2.plot(wave, (roall12[:,i]-riall12[:,i])/roall12[:,i], color=color[i])
plt.axis([ 265.9, 1035, -0.00333, 0.02540])
plt.xlabel("Wavelength [nm]", size=13)
plt.ylabel("Normalized ratio", size=13)
plt.tick_params(labelsize=12)
plt.grid()
plt.tight_layout()



===============================================================
                    Plot ratio_noise.png
===============================================================
rs= np.random.RandomState()
from datetime import datetime
print(str(datetime.now()))
rieccd = np.zeros((700000, 10))
roeccd = np.zeros((700000, 10))
for i in range(1,11):
    ri = np.load('e300_bigpsf_noslit/r'+str(i)+'/in.npy')
    ro = np.load('e300_bigpsf_noslit/r'+str(i)+'/out.npy')
    no = rs.poisson(np.mean(ro), size=(1000, 10240))
    ni = rs.poisson(np.mean(ri), size=(1000, 10240))

    ro, ri, wave, delta = the_thing((ro+no)*CCD, (ri+ni)*CCD)

    roeccd[:,i-1] = ro
    rieccd[:,i-1] = ri
ro_meaneccd = np.zeros((ro.shape[0]))
ri_meaneccd = np.zeros((ri.shape[0]))
for i in range(ri.shape[0]):
    ri_meaneccd[i] = np.mean(rieccd[i,:])
    ro_meaneccd[i] = np.mean(roeccd[i,:])
print(str(datetime.now()))

plt.figure(figsize=(9,5))
plt.plot(wave, (ro_meane-ri_meane)/ro_meane, color="tab:red", linewidth=3, label="Larger PSF")
#plt.plot(wave, (ro_meaneccd-ri_meaneccd)/ro_meaneccd, color="tab:green", linewidth=3, label="CCD included")
plt.plot(wave, (ro_mean_ccdnoise-ri_mean_ccdnoise)/ro_mean_ccdnoise, color="tab:blue", linewidth=3, label="Noise and CCD included")
plt.plot(in_spec[:,0], (in_spec[:,1]-in_spec2[:,1])/in_spec[:,1], color="tab:orange", linewidth=3, label="Input spectra ratio")
#plt.axis([ 265.9, 1035, -0.00333, 0.02540])
plt.xlabel("Wavelength [nm]", size=13)
plt.ylabel("Transit depth", size=13)
plt.grid()
plt.tick_params(labelsize=12)
plt.legend(fontsize=12)
plt.tight_layout()



===============================================================
                    Plot ratio_noise_err.png
===============================================================
err = np.zeros((700000))
for i in range(700000):
    err[i] = np.sqrt((1/10)* ((np.std(ro_ccdnoise[i,:])**2)/ro_mean_ccdnoise[i]**2 + (np.std(ri_ccdnoise[i,:])**2)/(ri_mean_ccdnoise[i]**2)))

div = (ro_mean_ccdnoise-ri_mean_ccdnoise)/ro_mean_ccdnoise

plt.figure(figsize=(9,5))
plt.plot(wave, (ro_meane-ri_meane)/ro_meane, color="tab:red", linewidth=3, label="Larger PSF")
#plt.plot(wave, (ro_meaneccd-ri_meaneccd)/ro_meaneccd, color="tab:green", linewidth=3, label="CCD included")
plt.plot(wave, (ro_mean_ccdnoise-ri_mean_ccdnoise)/ro_mean_ccdnoise, color="tab:blue", linewidth=3, label="Noise and CCD included")
plt.plot(in_spec[:,0], (in_spec[:,1]-in_spec2[:,1])/in_spec[:,1], color="tab:orange", linewidth=3, label="Input spectra ratio")
plt.fill_between(wave, div-err, div+err, color='k', alpha=0.2, label="Uncertainty")
#plt.plot(wave, div+err, 'k--', label="Uncertainty")
#plt.axis([ 265.9, 1035, -0.00333, 0.02540])
plt.xlabel("Wavelength [nm]", size=13)
plt.ylabel("Transit depth", size=13)
plt.grid()
plt.tick_params(labelsize=12)
plt.legend(fontsize=12)
plt.tight_layout()




"""

# print('Next image: \n')

# x_j, y_j = simfun.jitter(steps=int(exp*step), gain=0.4, amplitude_act=0.03, amplitude_sens=0.03) #New jitter, will have epx*step length
# jitter = np.stack((x_j, y_j), axis=-1)
# jitter =simfun.jitter_im(jitter[:,0], jitter[:,1], psf_size=(psf[:,:,0].shape[0], psf[:,:,0].shape[0]))

# image2, image2_wl=simfun.disperser(wl_endpoints=wl_ran, jit_img=jitter, psf_img=psf, pos=slitpos, image_size=img_size, 
#                                         dispersion=disper, eff=spec_eff2, mask_img=mask, steps=1, plot='n')
# # image_all = image + background

# img1 = simfun.bin_sum(image1, bin_size=sub_pixel)
# img2 = simfun.bin_sum(image2, bin_size=sub_pixel)
# read_wo = simfun.read_out(img_summed,1)
# plt.plot(read_wo)

# img2_summed = simfun.bin_sum(image2, bin_size=sub_pixel)
# read_w = simfun.read_out(img2_summed,1)
# plt.plot(read_w)
# plt.plot(np.arange(0, 750)+115, in_spec2[:,1]/in_spec[:,1])
# plt.plot(read_no/read_w)
# plt.axis((105, 770, 0.989, 0.98932))

'''
peak_ref, _ = find_peaks(in_spec[:,1]*(-1), threshold=4, prominence=200, distance=3)
plt.figure()
plt.plot(in_spec[:,0], in_spec[:,1])
plt.plot(in_spec[:,0][peak_ref], in_spec[:,1][peak_ref], "x")

rr = simfun.read_out(image)
peaks, _ = find_peaks(-1*rr, distance=20, prominence=500000)
plt.figure()
plt.plot(rr)
plt.plot(peaks, rr[peaks], 'x')

from scipy.stats import linregress
a, b, r, p, s = linregress(peaks, in_spec[:,0][peak_ref])
wave = a*np.arange(0, 10240)+b
plt.figure()
plt.plot(wave, rr)


#DET VAR HERTIL DU NÅEDE MED PLOTS. Der skal laves en del. Start med at loade de iamge-filer du lavede i går
# figsize=(6.2, 9.1)
fig = plt.figure()
gs1 = fig.add_gridspec(nrows=6, ncols=2)
ax1 = fig.add_subplot(gs1[0,0])
ax2 = fig.add_subplot(gs1[1,0])
ax3 = fig.add_subplot(gs1[2,0])
ax4 = fig.add_subplot(gs1[3,0])
ax5 = fig.add_subplot(gs1[4,0])
ax6 = fig.add_subplot(gs1[5,0])

ax1.imshow(image1)
ax1.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

ax2.imshow(image3)
ax2.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

ax3.imshow(image4)
ax3.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

ax4.imshow(image5)
ax4.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

ax5.imshow(image6)
ax5.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

ax6.imshow(image7)
ax6.tick_params(labelcolor='w', top=False, bottom=False, left=False, right=False)

ax11 = fig.add_subplot(gs1[0,1])
ax11.plot(image1[125,:], 'r-')
ax11.yaxis.tick_right()
ax11.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
ax11.tick_params(labelsize=12)
ax1.plot(np.linspace(0,999, 1000), np.ones(1000)*125, 'r-')


ax21 = fig.add_subplot(gs1[1,1])
ax21.plot(image3[125,:], 'r-')
ax21.yaxis.tick_right()
ax21.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
ax21.tick_params(labelsize=12)
ax2.plot(np.linspace(0,999, 1000), np.ones(1000)*125, 'r-')

ax31 = fig.add_subplot(gs1[2,1])
ax31.plot(image4[125,:], 'r-')
ax31.axis([291, 530, 88.57, 88.75])
ax31.yaxis.tick_right()
ax31.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
ax31.tick_params(labelsize=12)
ax3.plot(np.linspace(0,999, 1000), np.ones(1000)*125, 'r-')

ax41 = fig.add_subplot(gs1[3,1])
ax41.plot(image5[:,100], 'r-')
ax41.yaxis.tick_right()
ax41.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
ax41.tick_params(labelsize=12)
ax4.plot(np.ones(1000)*125, np.linspace(0,249, 1000), 'r-')

ax51 = fig.add_subplot(gs1[4,1])
ax51.plot(image6[:,100], 'r-')
ax51.yaxis.tick_right()
ax51.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
ax51.tick_params(labelsize=12)
ax5.plot(np.ones(1000)*125, np.linspace(0,249, 1000), 'r-')

ax61 = fig.add_subplot(gs1[5,1])
ax61.plot(image7[150, :], 'r-')
ax61.yaxis.tick_right()
ax61.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
ax61.tick_params(labelsize=12)
ax6.plot(np.linspace(0,999, 1000), np.ones(1000)*150, 'r-')

ax1.text(-65, 140, "A", size=18)
ax2.text(-65, 140, "B", size=18)
ax3.text(-65, 140, "C", size=18)
ax4.text(-65, 140, "D", size=18)
ax5.text(-65, 140, "E", size=18)
ax6.text(-65, 140, "F", size=18)


# ax7 = fig.add_subplot(gs1[6,0], autoscale_on=False)
# ax7.plot(simfun.read_out(image7))


fig = plt.figure()
ax = fig.add_subplot(111, frameon=False)
ax.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212, sharex=ax1)

temp = simfun.interp(np.loadtxt('optics.txt')[:,0], np.loadtxt('optics.txt')[:,1], wl_ran=wl_ran, kind="cubic", lowlim=wl_ran[0]-50, uplim=wl_ran[1]+50)
ax1.plot(np.loadtxt('optics.txt')[:,0],np.loadtxt('optics.txt')[:,1], label="Input")
ax1.plot(temp[:,0], temp[:,1], label="Interpolated")
ax1.tick_params(length=6, width=1, labelsize=12)
ax1.axis((300, 1337, 0.853404, 1))
ax1.grid()
ax1.legend(fontsize=12, loc=4)

temp = simfun.interp(np.loadtxt('QE.txt')[:,0], np.loadtxt('QE.txt')[:,1], wl_ran=[400,1250], kind="quadratic", lowlim=400, uplim=1250)
ax2.plot(np.loadtxt('QE.txt')[:,0],np.loadtxt('QE.txt')[:,1], label="Input")
ax2.plot(temp[:,0],temp[:,1], label="Interpolated")
ax2.tick_params(length=6, width=1, labelsize=12)
ax2.grid()
ax2.legend(fontsize=12, loc=1)
ax2.set_xlabel("Wavelength [nm]", size=12)
#ax.text(0,0.5,"Spectral throughput, $\eta$")
ax.text(-0.12, 0.5, 'Spectral throughput, $\eta$', va='center', rotation='vertical', size=12)
plt.setp(ax1.get_xticklabels(), visible=False)

# x_j, y_j = simfun.jitter(steps=int(3000), gain=0.7, amplitude_act=0.15, amplitude_sens=0.15)
plt.figure(figsize=(6.4, 6.4))
plt.plot(x_j, y_j)
plt.plot(np.linspace(-5,5,100), np.ones(100)*-5, 'k-')
plt.plot(np.linspace(-5,5,100), np.ones(100)*5, 'k-')
plt.plot(np.ones(100)*5, np.linspace(-5,5,100), 'k-')
plt.plot(np.ones(100)*-5, np.linspace(-5,5,100), 'k-')
plt.grid()
plt.xlabel('Sub-pixel', size=12)
plt.ylabel('Sub-pixel', size=12)
plt.axis((-12,12,-12,12))
'''