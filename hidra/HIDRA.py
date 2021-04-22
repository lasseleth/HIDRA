import numpy as np
import sys

def same_dist_elems(arr):
    """
    Smart little script to check if indices are equidistant. 
    Found at https://stackoverflow.com/questions/58741961/how-to-check-if-consecutive-elements-of-array-are-evenly-spaced
    
    Parameters
    ----------
    arr : array_like
        Input array
        
    Returns
    -------
    bool
        boolean value, True if array is equidistantly spaced, False otherwise
    """
    diff = arr[1] - arr[0]
    for x in range(1, len(arr) - 1):
        if arr[x + 1] - arr[x] != diff:
            return False
    return True

def progressbar(it, prefix="", size=60, file=sys.stdout):
    """
    Function to generate a progress bar. Does not work ideally... Found on stackexchange
    """
    count = len(it)
    def show(j):
        x = int(size*j/count)
        file.write("%s[%s%s] %i/%i\r" % (prefix, "#"*x, "."*(size-x), j, count))
        file.flush()        
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    file.write("\n")
    file.flush()

def mult(*args):    
    # Multiply elements one by one 
    result = 1
    for x in args: 
         result = result * x  
    return result

def interp(x, y, wl_ran=(300, 1200), delta_lambda=1, kind='cubic', lowlim=400, uplim=1100):
    """
    This function interpolates values between given input table values.
    
    Parameters
    ----------
    x : array_like
        Input array of x values, eg. wavelength
        Ex: np.array([100, 217, 350])
    y : array_like
        Input array of y values, eg. quantum efficieny, or mirror reflectance. 
        Ex: np.array([0.1, 0.7, 0.85]) 
    wl_ran : tuple
        wavelength span. Entries must be integers
    delta_lambda : float
        wavelength resolution, in nm.
    kind : string
        type of interpolation. Valid options are 'linear', 'quadratic' and 'cubic'. 
    lowlim : float
        lower wavelength limit. Below this value, throughput will be set to 0
    uplim : float
        upper wavelength limit. Above this value, thoughput will be set to 0
        
    Returns
    -------
    interpolated : array_like
        Interpolated values between wl_start and wl_stop, with sharp cutoff beyond the specified limits. 
        
    Notes
    -----    
    Check interp1d for more options.
    
    """
    from scipy.interpolate import interp1d #Load neccessary package
    import numpy as np
    f = interp1d(x, y, kind=kind, fill_value="extrapolate") #interpolates, and extrapolates if the given table does not cover the wavelength range
    # xnew = np.linspace(wl_ran[0], wl_ran[1], num=int((wl_ran[1]-wl_ran[0])/delta_lambda), endpoint=True) #Generates new x-values
    xnew = np.arange(wl_ran[0], wl_ran[1], delta_lambda)
    interp = f(xnew) #"Raw" interpolation
    interpol= np.asarray([i if i>0 else 0 for i in interp]) #recast as numpy array for easier handling, and throws away values below 0
    interpolated = np.stack((xnew,interpol), axis=-1) #Combine new x-values and interpolated
    
    # To remove values below lower limit
    for i in range(interpolated.shape[0]):
        if interpolated[i,0]<lowlim:
            interpolated[i,1]=0
        if interpolated[i,0] > lowlim:
            break
    
    #To remove values above upper limit
    for i in reversed(range(interpolated.shape[0])): #Start from top and goes down
        if interpolated[i,0]>uplim:
            interpolated[i,1]=0
        if interpolated[i,0] < uplim:
            break
    
    return interpolated

def loadfunc(*args, wls, kind='cubic'):
    result = 1
    for x in args: #takes several input arrays
        loaded = np.loadtxt(x)
        if not loaded.shape[0] == (wls[1]-wls[0]): #if input is not of the correct length, this will interpolate
            temp = interp(loaded[:,0], loaded[:,1], wl_ran=wls, kind=kind, lowlim=wls[0]-50, uplim=wls[1]+50)[:,1]
        else:
            temp = loaded
        result = result * temp
    return result

def CCD_maker(CCD_size, subpix=10, var=0.05, var2=0.05, grid_loss=0.6, smooth=5):
    """
    This function creates a CCD composed of subpixels, with a separating grid between all full pixels. 
    The grid will have some loss.
    
    Parameters
    ----------
    CCD_size : array_like
        Input size of CCD (in full pixels). 
        Ex: (10, 10)
    subpix : int
        Number of subpixels in each full pixel
    var : float
        Variation of noise, (in the gaussian noise)
    var2 : float
        Variation, relative variation from 0
    grid_loss: float
        Loss in the grid. 1 = everything gets through, 0 = nothing gets through
    smooth : float
        Smoothness factor, previously called "stepsize". Is the ammount of subpixel to correlate to during that phase. Must be larger than 1.
        
    Returns
    -------
    CCD : ndarray
        Output array of the CCD with the specified subpixel ammount, and size.
        
    Notes
    -----
    Once used, remember to save the created CCD, as to not run the script again. It can take quite a while to make big arrays.

    Exaples
    -------
    >>> new_CCD = CCD_maker((240, 240), 10, 0.3, 0.7, 5)
    array([[0.59858663, 0.59919131, 0.59980866, ..., 0.59164421, 0.59224492,
        0.59108706],
       ...,
       [0.63641557, 0.88710319, 0.60372464, ..., 0.91472067, 0.65503371,
        0.96646196]])    
    """
    import numpy as np
    import sys
    gridsize = subpix #"size of pixels" in # of subpixels
    x_size = CCD_size[0]*gridsize
    y_size = CCD_size[1]*gridsize # number of subpixels 
    S = smooth #stepsize "smoothness", previously =5
    
    CCD = np.random.normal(1-var, var, [x_size, y_size])#*var #Noise matrix 
    # noise = np.random.standard_normal((x_size, y_size))*var #Noise matrix 
    # CCD = np.ones((x_size,y_size)) #"Clean" matrix
    # CCD = CCD-noise #Subtracts noise from "clean" CCD matrix
    CCD2=np.zeros(CCD.shape)
    
    #Correlate the subpixels
    N = 3 # number of times to correlate
    for t in np.arange(0,N):
        for i in np.arange(0,x_size):
            for j in np.arange(0,y_size):
                bit = CCD[i:i+S, j:j+S-1] #cuts out a bit to treat
                CCD2[i, j] = np.sum(np.sum(bit)/np.size(bit)) #correlates surrounding subpixels
        sys.stdout.write('.'); sys.stdout.flush(); #"Progress bar", just for visuals   
    #Introduces grid, to mimic the actual pixels - seperate the subpixels by a grid with a slight loss defined by grid_loss variable. 
    grid = np.ones((CCD.shape[0], CCD.shape[1])) #Set up grid
    grid[0::gridsize,:]=grid_loss #Sets gridloss for every 'gridsize' row (10)
    grid[:,0::gridsize]=grid_loss #sets gridloss for every 'gridsize' coloumn (10)
    #to see a visualization of this, use the variable explorer - type: %varexp --imshow grid
    # noise2 = np.random.standard_normal((x_size, y_size))*var2 
    noise2 = np.random.normal(0, var2, [x_size, y_size])#*var2 
    # CCD2 = CCD2+noise2+1
    CCD2 = CCD2-noise2
    
    CCD2 = CCD2/np.mean(CCD2)
    # CCD2 = CCD2/np.mean(CCD2)
    CCD = CCD2*grid #overlays the grid on the CCD
    # CCD = CCD/np.max(CCD)
    
    return CCD

def psf_maker(file_name, wl_endpoints=(350, 1100), f=1, size=101, res=(100, 100)):
    """
    Creates a new file containing the full-color PSF, without interpolating as with psf_maker

    Parameters
    ----------
    file_name : str
        Desired name of the file.
    wl_endpoints : tuple, optional
        Two values that mark the first and last colors. The default is (350, 1100).
    f : float
        factor to multiply in the sigma values
    size : int, optional
        Size of the PSF. The default is 101.
    res : tuple, optional
        Resolutionof the meshgrid used in the 2D Gaussian. Will affect the size of the PSF inversly: Larger values mean smaller PSF. Just a tweakable parameter. The default is (100, 100).

    Returns
    -------
    .npy and .hdf5 files containing the PSF
    """
    
    import os
    import numpy as np
    import h5py
    
    path = os.getcwd() #Get current working directory
    file_path = path + "/" + file_name +".hdf5" #Set up path to save file later
    ''' 
    numColors = int( (wl_endpoints[1]-wl_endpoints[0])/step) # Number of colors
    x_size = size[0]
    y_size = size[1] #Extracts from the size input
    
    z = np.float128(np.zeros((res, res, numColors))) # Setup empty array for PSF-slices
    
    x = np.float128(np.linspace(-x_size, x_size, res)) #Preparation for meshgrid
    y = np.float128(np.linspace(-y_size, y_size, res))
    xx, yy = np.meshgrid(x, y) #define meshgrid
    
    for i in range(wl_endpoints[0], wl_endpoints[1], step): # for-loop to create one psf for each color
        sigma_x = np.float128(np.log(i)+0.5*i/100) # Used in the 2D Gaussian
        sigma_y = np.float128(np.log(i)+0.5*i/100)
    
        # 2D Gaussian function, that takes sigma_x and _y as input variables. Also takes in the meshgrid xx and yy
    
        zz = (1/(2*np.pi*sigma_x*sigma_y) * np.exp(-((xx)**2/(2*sigma_x**2) 
             + (yy)**2/(2*sigma_y**2))))    
    
        zz = zz/np.sum(zz) # Normalizes, so the total value (the sum of the array) =1
        z[:,:,i-350] = zz # put psf-"slice" into larger 3D array    
    
    '''
    step=1
    numColors = int( (wl_endpoints[1]-wl_endpoints[0])/step) # Number of colors
    x_size = res[0]
    y_size = res[1] #Extracts from the size input
    
    z = np.zeros((size, size, numColors)) # Setup empty array for PSF-slices
    
    x = np.linspace(-x_size, x_size, size) #Preparation for meshgrid
    y = np.linspace(-y_size, y_size, size)
    xx, yy = np.meshgrid(x, y) #define meshgrid
    
    for i in range(wl_endpoints[0], wl_endpoints[1], step): # for-loop to create one psf for each color
        # sigma_x = np.log(i)+f*i/100 # Used in the 2D Gaussian, old one
        # sigma_y = np.log(i)+f*i/100
        
        sigma_x = f*0.014285714285714285 * i + 20.714285714285715 # emperically determined slope, linear increase 
        sigma_y = f*0.014285714285714285 * i + 20.714285714285715 
        
        # 2D Gaussian function, that takes sigma_x and _y as input variables. Also takes in the meshgrid xx and yy
    
        zz = (1/(2*np.pi*sigma_x*sigma_y) * np.exp(-((xx)**2/(2*sigma_x**2) 
             + (yy)**2/(2*sigma_y**2))))    
        
        zz = zz/np.sum(zz) # Normalizes, so the total value (the sum of the array) =1
        z[:,:,i-wl_endpoints[0]] = zz # put psf-"slice" into larger 3D array
    
    if os.path.exists(file_path) == True: #If file already exists, it will be deleted
        os.remove(file_path)
        
    # Saving the psf as a hdf5 file in order to store the large file, using h5py
    psf_file = h5py.File(file_path, "a")
    psf_file.create_dataset('psf', data=z, dtype='f') # Place dataset in the .hdf5 file

    np.save(file_name + "_raw.npy", z) #Save as .npy binary file

    return print("New PSF done, saved as", file_name, ".npy")

def psf_interp(input_psf_images, input_psf_wl, wl_endpoints=(350, 1100), delta_lambda=1):    
    import sys
    from scipy.interpolate import interp1d
    print('\nInterpolating missing wavelengths in PSF... \n')
    ran = range(wl_endpoints[0], wl_endpoints[1], delta_lambda) #set for-loop range
    res = input_psf_images.shape[0] # Width of the input psf, so the created psf will have the same size

    psf = np.zeros((input_psf_images.shape[0], input_psf_images.shape[1], wl_endpoints[1]-wl_endpoints[0])) #Creates empty array for the new psf
    for i in range(res):
        for j in range(res):
            f_test = interp1d(input_psf_wl, input_psf_images[i,j,:], kind='quadratic', fill_value="extrapolate") #sets up interpolation function
            psf[i,j,:] = f_test(ran) # interpolates at the wavelengths specified in the range
        sys.stdout.write('.'); sys.stdout.flush(); #"Progress bar", just for visuals
    print(' ')
    print('Interpolation done')
    print(' ')    
    return psf

def func_jitter (entries, gain, dt):
    """
    Generates two jitter arrays, in x- and y. 
    
    Parameters
    ----------
    entries : int
        Number of entries in the desired jitter arrays
    gain : float
        Gain of the ADCS. 
    dt : int
        Time delay

    Returns
    -------
    x, y : array-like
        Jitter in x- and y-directions
    """
    x = np.zeros((entries+dt)) #allocates for arrays
    y = np.zeros((entries+dt))
    
    for i in range(entries+dt-1): #set up for loop
        x[i+1] = x[i]+np.random.normal()-gain*x[i-dt] #next entry will be previous, plus a Gaussian number, 
        y[i+1] = y[i]+np.random.normal()-gain*y[i-dt] # and the correction is subtracted from the i-dt'th entry
    x = x[dt-1:-1] #Cut off the initial dt entries. 
    y = y[dt-1:-1]
    return x, y

def func_slit(slit_size=[10,100], pos=[499, 499], image_size=[1000,1000]):
    """ Creates a slit "mask" to overlay images.
    
    Parameters
    ----------
    slit_size : array_like, int
        Size of slit: should be two numbers, width and height.
    pos : array_like, int
        Position of the slit, measured in subpixels.
    img_size : array_like, int
        Size of mask. Should be identical to size of the image upon which the mask is overlaid.
        
    Returns
    -------
    mask : array_like
        Mask is zero everywhere except in the slit, where the value is 1.
    """
    width = slit_size[0] #Loads in size of slit
    height = slit_size[1]
    
    x_low = pos[0] - width #Finds boundaries
    x_up = pos[0] + width
    y_low = pos[1] - height
    y_up = pos[1] + height
    
    mask = np.zeros(image_size) #Creates empty mask
    
    mask[y_low:y_up, x_low:x_up] = mask[y_low:y_up, x_low:x_up]+1 #Fills in the slit, so that only the slit has any throughput
    
    return mask

def mag(mag_star, mag_ref=0):
    """Calculates the brightness difference based on magnitudes
    Parameters
    ----------
    mag_star : float
        Magnitude of input star
    mag_ref : float
        magnitude of reference star"""
    return 10**(0.4*((mag_ref)-(mag_star)))



# def jitter_im(x, y, psf_size):
#     ''' Creates a jitter "image" - a matrix of the same dimensions (x & y) as the psf, used in the folding function
#     NOTE: Will round of the position of the jitter to nearest subpixel!

#     Parameters
#     ----------
#     x : array
#         Input jitter x-coord.
#     y : array
#         Input jitter y-coord.
#     psf_size : int, two values
#         Size of the psf.

#     Returns
#     -------
#     jitter : array
#         Jitter image, where each point where the jitter "stops" has a +1 value. All other points are zero.
#     '''
#     jitter=np.zeros(psf_size) # Setup image
#     # jitter2=np.zeros(psf_size)
#     for i in range(len(x)):
#         jitter[(x[i]+(psf_size[0]/2)).astype(int), (y[i]+(psf_size[1]/2)).astype(int)]= jitter[(x[i]+(psf_size[0]/2)).astype(int), (y[i]+(psf_size[1]/2)).astype(int)]+1
#         # jitter2[x[i].astype(int)+int(np.floor(psf_size[0]/2)), y[i].astype(int)+int(np.floor(psf_size[1]/2))]= jitter[x[i].astype(int)+int(np.floor(psf_size[0]/2)), y[i].astype(int)+int(np.floor(psf_size[1]/2))]+1 # Create jitter "image". +1 to every point where the jitter "hits"
#     return jitter#, jitter2

def jitter_im(x, y, psf_size):
    ''' Creates a jitter "image" - a matrix of the same dimensions (x & y) as the psf, used in the folding function
    NOTE: Will round of the position of the jitter to nearest subpixel!

    Parameters
    ----------
    x : array
        Input jitter x-coord.
    y : array
        Input jitter y-coord.
    psf_size : int, two values
        Size of the psf.

    Returns
    -------
    jitter : array
        Jitter image, where each point where the jitter "stops" has a +1 value. All other points are zero.
    '''
    jitter=np.zeros(psf_size) # Setup image
    # jitter2=np.zeros(psf_size)
    for i in range(len(x)):
        rang1 = (x[i]+(psf_size[0]/2)).astype(int) 
        rang2 = (y[i]+(psf_size[1]/2)).astype(int)
        # print(rang1, rang2)
        jitter[rang1, rang2] = jitter[rang1, rang2]+1
        # Create jitter "image". +1 to every point where the jitter "hits"
    return jitter#, jitter2


def folding(psf_image, jitter_image, mode='same', boundary='fill'):
    #Clutter function, might as well just use signal.convolve2d 
    from scipy import signal
    folded=signal.convolve2d(psf_image, jitter_image, mode=mode, boundary=boundary) #convolves the psf slice and the jitter image
    return folded

#The correct disperser::::
def disperser(wl_endpoints, jit_img, psf_ends, pos, image_size, dispersion, eff, 
              mask_img, steps=1, secondary_source='n', plot='n'):
    import sys
    from scipy import signal
    from astropy.convolution import AiryDisk2DKernel
    x_pos=pos[0] 
    y_pos=pos[1] #load in position of "zeroth order"
    im_disp = np.zeros((image_size[0],image_size[1])) # empty image
    im_disp_lambda = np.zeros((image_size[0],image_size[1])) 
    x_dispersion = dispersion[0] #load in dispersions
    y_dispersion = dispersion[1]
    numColors = int( (wl_endpoints[1]-wl_endpoints[0])) #total number of colours to iterate
    # print("Number of colors to iterate: " + str(numColors))
    # print(' ')
    if plot=='y': #this part is not useful atm
        import matplotlib.pyplot as plt
        plt.figure()
        from matplotlib.colors import LinearSegmentedColormap
        N = 256 #8-bit value, to fix colours
        colspec = plt.cm.get_cmap('Spectral') #Fetches colourmap to use later
        vals = np.ones((N,4)) #Setup for colormap
    temp = np.linspace(psf_ends[0], psf_ends[1], numColors)
    for i in range(0, numColors, steps):
    # for i in range(0,101, steps):
        im = np.zeros((image_size[0],image_size[1])) #create temp. image
        
        psf = AiryDisk2DKernel(temp[i]*0.1, x_size=jit_img.shape[0], y_size=jit_img.shape[0]).array #PSF for this colour
        
        if secondary_source == 'y': #To account for the secondary light source perhaps not being fully within the psf
            # fold = folding(psf_img[:,:,i], jit_img)
            fold = signal.convolve2d(psf[:,:,i], jit_img, mode='same', boundary='fill') #fold psf and jitter
            fold = fold[0:jit_img.shape[1], 0:jit_img.shape[0]] #cut down to regular shape
        else:
            fold = signal.convolve2d(psf[:,:], jit_img, mode='same', boundary='fill') #fold as usual, if no sec. sources
        # fold=fold/np.sum(fold)
        
        foo = int(psf.shape[0]/2)
        # im[0+x_pos-foo:len(jitter)+x_pos-foo, 0+y_pos-foo:len(jitter)+y_pos-foo] = im[0+x_pos-foo:len(jitter)+x_pos-foo, 0+y_pos-foo:len(jitter)+y_pos-foo] + fold*magni
        im[0+y_pos-foo:len(fold)+y_pos-foo, 0+x_pos-foo:len(fold)+x_pos-foo] = fold #im[0+y_pos-foo:len(fold)+y_pos-foo, 0+x_pos-foo:len(fold)+x_pos-foo] + fold#*magni
        immask = im*mask_img #mask is "overlaid" by multiplying
        
        roll_x = np.roll(immask,  int(np.modf(x_dispersion[i])[1]), axis=1) #move/disperse the light 
        roll_y = np.roll(roll_x,  int(np.modf(y_dispersion[i])[1]), axis=0) #also in the y-direction

        dx = abs(np.modf(x_dispersion[i])[0]) #residual amount (decimal amounts are shifted to the next sub-pixel)
        dy = abs(np.modf(y_dispersion[i])[0])
        
        foob = roll_y*(eff[i]*(1-dx)*(1-dy)) #multiply by efficiency 
        im_disp = im_disp + foob  # Add the rolled image to the final, and multiply by the "effectivity"
        
        
        roll_dx = np.roll(roll_y, 1, axis=1) # Roll the residual to the next subpixel
        eff_dx = eff[i] * dx * (1-dy) # effectivity of the x-residual
        
        roll_dy = np.roll(roll_y, 1, axis=0) # Roll the residual to the next subpixel, y-wise
        eff_dy = eff[i] * dy * (1-dx) # y-residual eff.
        
        roll_dxy = np.roll(roll_dx, 1, axis=0) # roll the image one step in both x- and y-wise.
        eff_dxy = eff[i]* dx * dy   #and eff.
        
        baar = roll_dx*eff_dx + roll_dy*eff_dy + roll_dxy*eff_dxy
        
        im_disp = im_disp + baar #add all residuals and multiply by their respective effectivities.
        
        im_disp_lambda = im_disp_lambda+((foob+baar)*(i+wl_endpoints[0])) #fill in im_disp, and multiply by wavelength i
        # im_disp_lambda = im_disp_lambda+(i+wl_endpoints[0]) #fill in im_disp, and multiply by wavelength i
        # sys.stdout.write('/'); sys.stdout.flush(); #"Progress bar", just for visuals    
        
        ##### Plotting #####
        if plot == 'y':
            vals[:, 0] = np.linspace(0, colspec(1-i/750)[0], N) #Making new colourmap values
            vals[:, 1] = np.linspace(0, colspec(1-i/750)[1], N) #the /750 is to normalize the colormap, so values fall between 0 and 1
            vals[:, 2] = np.linspace(0, colspec(1-i/750)[2], N)
            vals[:, 3] = np.linspace(0, 1, N) #alpha, for making the cmap transparent
            newcmp = LinearSegmentedColormap.from_list(name='Spectral', colors=vals) #Creates new cmp, based on vals    
            plt.imshow(roll_y, cmap=newcmp) # Show array   
    
    if plot=='y':
        plt.title('Color dispersion of sample spectrum', size=18)
        plt.xlabel('Sub-pixel', size=13)
        plt.ylabel('Sub-pixel', size=13)
    return im_disp, im_disp_lambda






def ccd_interp(inCCD, wls, img, img_wl):
    """
    Interpolator used to find the subpixel sensitivity for all wavelengths (not just the ones created by ccd_maker)

    Parameters
    ----------
    inCCD : array
        Input CCD array, can be made using ccd_maker.
    wls : array
        Corresponding wavelengths. Must have the same size as the depth of inCCD
    img : array
        Input image, from disperser2.
    img_wl : array
        Input image wavelengths.

    Returns
    -------
    new_img : array
        Image "multiplied" by the CCD, using the interpolated sensitivities for each subpixel.
    """
    import sys
    from scipy.interpolate import interp1d
    
    if not wls.shape[0] is inCCD.shape[2]:
        raise TypeError("Wavelength array and input CCD depth not same size")
    if not inCCD.shape[0:2] == img.shape[0:2] == img_wl.shape:
        raise TypeError("CCD and image not same size")    
    
    new_img = np.zeros((img.shape[0], img.shape[1]))
    for i in range(0, inCCD.shape[0]):
        for j in range(0, inCCD.shape[1]):
            interp = interp1d(wls, inCCD[i,j,:], kind="slinear", fill_value="extrapolate")
            new_img[i,j] = img[i,j]*interp(img_wl[i,j])
        sys.stdout.write('.'); sys.stdout.flush();
    return new_img


def read_out(dispersed):
    '''
    Will sum up the "photons" in the y-direction of the input dispersed image.
    
    Parameters
    ----------
    dispersed : array, 2 dimensional
        Dispersed image-array.

    Returns
    -------
    counts : array
        Array of counts in the y-direction.

    '''
    counts = np.array(())
    for i in range(dispersed.shape[1]):
        counts = np.append(counts, np.sum(dispersed[:,i]))
    return counts
def read_outx(dispersed):
    '''
    Will sum up the "photons" in the X-direction of the input dispersed image.
    '''
    counts = np.array(())
    for i in range(dispersed.shape[0]):
        counts = np.append(counts, np.sum(dispersed[i,:]))
    return counts

def bin_sum(inp, bin_size): 
    """
    Returns a binned version of inp, with each bin being bin_size in each dimension. The bins are summed up.

    Parameters
    ----------
    inp : array_like
        Input array. Must be 2D.
    bin_size : int
        Bin size. Division of input shape and bin_size should be a whole number, i.e. no 8.333 etc.

    Returns
    -------
    binned : array
        Array of inp.shape/bin_size in shape, with the bins summed up.

    """
    # Check if bin_size is whole divisor of inp.shape
    if not np.modf(inp.shape[0]/bin_size)[0] == 0 == np.modf(inp.shape[1]/bin_size)[0]:
        raise TypeError("Input shape and bin size divided must be a whole number. (mod = 0)")
    
    temp = np.zeros((inp.shape[0], int(inp.shape[1]/bin_size) )) #Create empty matrix for first step
    summed = np.zeros((int(inp.shape[0]/bin_size), int(inp.shape[1]/bin_size) )) #Empty matrix for second step

    for x in range(0, inp.shape[1], bin_size): #Range for 1st
        j = range(0+x, bin_size+x) #Bin range. ex. 20-30 if bin_size is 10
        for i in range(0, inp.shape[0]): # over all columns 
            temp[i, int(j[0]/bin_size)]= sum(inp[i,j]) #sum, and add to temp
    
    for x in range(0, inp.shape[0], bin_size): #2nd step, repeat 1st step, but for rows
        i = range(0+x, bin_size+x) #row bin-range. 
        for j in range(0, summed.shape[1]): 
            summed[int(i[0]/bin_size), j]= sum(temp[i,j]) #sum and add to result-matrix
            
    return summed

def noise(size, image, RON=5):
    noise = np.zeros((size[0], size[1]))
    for i in range(size[0]):
        for j in range(size[1]):
            noise[i,j] = (np.sqrt(image[i,j])+RON)*np.random.normal(0,1)
    return noise

def convert_plate_pix(plate_scale, pix_size):
    """
    Plate scale is calculated with the equation:
    P = 206265 / (D*f/#)
    206265 is the amount of arcsecs in a radian.
    D is the diameter of the telescope
    f/# is the f-number: Focal length/Diameter
    ( http://www-supernova.lbl.gov/~sed/telescope/obsguide/platescale.html )
    
    For a telescope of 20 cm, and focal length of 50 cm, the plate scale is 412.53 arcsec/mm
    Parameters
    ----------
    plate_scale : float
        Must be in arcsec/mm.
    pix_size : float
        Must be in mm/pixel.

    Returns
    -------
    convert_factor : float
        How large a sky area a single pixel width covers.
    """
    convert_factor = plate_scale * pix_size # [arcsec per pix] = [arcsec/mm] * [mm/pix] 
    return convert_factor

def convert_slit(unit, size, convert_factor):
    if not type(unit) == str:
        raise TypeError("unit must be a string")
    if not ((unit == 'pix') or (unit == 'ang')):
        raise TypeError("unit must be either 'ang' or 'pix'")
    
    if unit == 'ang':
        slit_size = np.divide(size, convert_factor)
    if unit == 'pix':
        slit_size = size
        
    return slit_size

def setup(input_file):
    import warnings
    in_spec = np.loadtxt(input_file.in_spec) 	#Input science spectrum
    in_spec2 = np.loadtxt(input_file.in_spec2) 	#Input science spectrum
    
    col_area = input_file.col_area	    	# Collecting area
    
    sub_pixel = input_file.sub_pixel       # Amount of sub-pixels per full pixel
    
    img_size = input_file.img_size 		#Size of the CCD, in pixels
    
    pl_scale = input_file.pl_scale 		# Plate scale
    
    pix_size = input_file.pix_size 		# Pixel size
    
    bg_spec = np.loadtxt(input_file.bg_spec) 	# Background spectrum, i.e. zodiacal light
    
    exp = input_file.exp 			       	# Exposure time
    
    wl_ran = input_file.wl_ran 			# Wavelength range
    
    eta_in = input_file.eta_in 	# Spectral troughput of the entire system. Requires at minimum the CCD QE
    
    slit = input_file.slit 	           	# Slit size. Unit first, then width and height
    
    #psf = np.load(input_file.psf) 		        	# Point Spread Function of the optics etc.
    
    #psf_col = input_file.psf_col
    #psf_col = np.arange(300, 1000)
    
    disper = np.load(input_file.disper) 	#Dispersion of the spectrograph
    
    ####### Optionals ########
    if not input_file.jitter:
        jitter = ''
    else:
        jitter = np.load(input_file.jitter) 			     	# Spacecraft jitter   
    
    
    step = input_file.step 				# Step size. Only needed if jitter is left empty
    
    in_CCD = np.load(input_file.in_CCD) 	   	# Input CCD imperfections. Sub-pixel variations
    
    CCD_col = input_file.CCD_col           # CCD colours, respective to each slice in in_CCD
    
    
    
    img_size[0] = img_size[0]*sub_pixel
    img_size[1] = img_size[1]*sub_pixel
    pl_arc_mm = convert_plate_pix(pl_scale, pix_size=pix_size)
    
    disper[0] = disper[0]*sub_pixel
    # disper[1] = disper[1]*sub_pixel
    
    
    
    if not ((type(wl_ran[0]) == int) or (type(wl_ran[1]) == int)):
        raise TypeError("wl_ran must be a tuple with two integers")
    span = wl_ran[1]-wl_ran[0]
    
    foo = 1
    args = eta_in
    for x in args: #takes several input arrays
        loaded = np.loadtxt(x)
        if not loaded.shape[0] == span: #if input is not of the correct length, this will interpolate
            temp = interp(loaded[:,0], loaded[:,1], wl_ran=wl_ran, kind='cubic', lowlim=wl_ran[0]-50, uplim=wl_ran[1]+50)[:,1]
        else:
            temp = loaded
        foo = foo * temp
    
    eta_in = foo
    del foo, args, temp, loaded
    
    
    
    #Handling the input spectrum and SEC/TEC
    if not in_spec.shape[0] == span:
        in_spec = interp(x=in_spec[:,0], y=in_spec[:,1], wl_ran=wl_ran, kind='cubic', lowlim=wl_ran[0]-50, uplim=wl_ran[1]+50)
    if not eta_in.shape[0] == span:
        raise TypeError("eta_in must cover the range of wavelengths: " + str(span) + " entries, from " + str(wl_ran[0]) +" to " +str(wl_ran[1]))
    spec_eff = in_spec[:,1] * col_area * eta_in
    
    if not in_spec2.shape[0] == span:
        in_spec2 = interp(x=in_spec2[:,0], y=in_spec2[:,1], wl_ran=wl_ran, kind='cubic', lowlim=wl_ran[0]-50, uplim=wl_ran[1]+50)
    if not eta_in.shape[0] == span:
        raise TypeError("eta_in must cover the range of wavelengths: " + str(span) + " entries, from " + str(wl_ran[0]) +" to " +str(wl_ran[1]))
    spec_eff2 = in_spec2[:,1] * col_area * eta_in
    
    
    #Slit is created here
    slit_size = convert_slit(unit = slit[0], size = slit[1:3], convert_factor = pl_arc_mm) #Convert slit size to pixels
    slit_size[0] = slit_size[0]*sub_pixel #Convert to subpixels
    slit_size[1] = slit_size[1]*sub_pixel
    slitpos = [150, 249] #Slit position on the sub-pixel CCD image. Arbitrary position..
    mask = func_slit(slit_size = np.floor(slit_size).astype(int), pos=slitpos, image_size=img_size) #Generate mask used to overlay before actual dispersion later.
    
    
    #Background spectrum gets handled here. A background image of exp = 1s will be created, and can be scaled and overlaid on the final image
    # new_bg = input("Do you wish to generate a new background? (y/n): ")
    new_bg = "n"
    if new_bg == 'y':
        if not bg_spec.shape[0] == span: #interpolate if values are missing
            bg_spec = interp(x=bg_spec[:,0], y=bg_spec[:,1], wl_ran=wl_ran, kind='cubic', lowlim=wl_ran[0]-50, uplim=wl_ran[1]+50)
            print("\nInterpolated missing values in background spectrum")
        detector_area = (pl_arc_mm*img_size[0]/sub_pixel)*(pl_arc_mm*img_size[1]/sub_pixel) #Collecting area of the detector measured in arcsec^2
        bg_spec = bg_spec*detector_area #Multiply by detector area
        bg_psf = np.ones((101, 101, wl_ran[1]-wl_ran[0]))
        x_j, y_j = func_jitter(entries=(exp*step), gain=0.15, dt=5) #This jitter will be a single point at the center of the jitter image
        bg_jit = jitter_im(x= x_j, y= y_j, psf_size=(bg_psf[:,:,0].shape[0], bg_psf[:,:,0].shape[0]) ) #Creating jitter "image"
        background, background_wl = disperser(wl_endpoints=wl_ran, jit_img=bg_jit, psf_img=bg_psf, pos=slitpos, image_size=img_size, dispersion=disper, eff = bg_spec[:,1], mask_img=mask, steps=1, plot='n' )
        np.save('background.npy', background) #saving the background image for later use. 
        del x_j, y_j, bg_jit, background_wl, bg_spec #getting rid of unnecessary variables
    else:
        background = np.load('../sample_values/background.npy')
    
    
    try: #If jitter is not defined, new jitter will be generated
        jitter
    except NameError:
        try:
            step
        except NameError:
            raise NameError("Either jitter or step must be specified")
        x_j, y_j = func_jitter(entries=(exp*step), gain=0.15, dt=5) 
        # x_j, y_j = simfun.jitter(entries=(exp*step), gain=0.02, dt=10)
        jitter = np.stack((x_j, y_j), axis=-1)
        spec_eff = spec_eff/step #If the generated jitter is used, the spectrum must be in step size, not seconds
        spec_eff2 = spec_eff2/step
    with warnings.catch_warnings(): #This is to suppress the potential "FutureWarning" error message. Comparing np-array to str etc. Might cause errors down the line?
        warnings.simplefilter(action='ignore', category=FutureWarning)
        if jitter == '': #if jitter is an empty str, it will also be generated.
            if step == '': #step must be specified
                raise TypeError("If jitter is unspecified, step must be explicitly specified")
            x_j, y_j = func_jitter(entries=(exp*step), gain=0.15, dt=5) #New jitter, will have epx*step length
            # x_j, y_j = simfun.jitter(entries=(exp*step), gain=0.02, dt=10)
            jitter = np.stack((x_j, y_j), axis=-1)
            spec_eff = spec_eff/step #If the generated jitter is used, the spectrum must be in step size, not seconds
            spec_eff2 = spec_eff2/step
    jitter = jitter_im(x= jitter[:,0], y= jitter[:,1], psf_size=(101, 101) )

    return spec_eff, spec_eff2, jitter, x_j, y_j, img_size, sub_pixel, pl_arc_mm, disper, mask, slitpos, background


def int_r(r1, rang):
    """
    Interpolator
    """
    from scipy.interpolate import interp1d
    x= np.arange(len(r1))
    xnew = np.arange(0, len(r1), 0.001)
    f1 = interp1d(x, r1, kind=3, fill_value="extrapolate")
    # f2 = interp1d(x, r2, kind=3, fill_value="extrapolate")
    
    r1_int = f1(xnew)
    # r2_int = f2(xnew)
    return r1_int

# def int_r(r1, r2, rang):
#     """
#     Interpolator
#     """
#     from scipy.interpolate import interp1d
#     x= np.arange(len(r1))
#     xnew = np.arange(0, len(r1), 0.001)
#     f1 = interp1d(x, r1, kind=3, fill_value="extrapolate")
#     f2 = interp1d(x, r2, kind=3, fill_value="extrapolate")
    
#     r1_int = f1(xnew)
#     r2_int = f2(xnew)
#     return r1_int, r2_int

def noise1d(x, RON=5):
    noise = np.zeros((x.shape))
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            noise[i,j] = (np.sqrt(x[i,j])+RON)*np.random.normal(0,1)
    return noise

def noise_test(x):
    noise = np.sqrt(abs(x))*np.random.normal(0, 1)
    return noise

def sinusoidal(size, frequency, amplitude, phase):
    """
    Function to generate sinusoidal jitter
    
    Parameters
    ----------
    size : int
        Size of the desired output array.
    frequency : list
        List with frequencies
    amplitude : list
        List of amplitudes. Must be the same size as the frequency list
    phase : float
        Phase of the sinusoidal.

    Returns
    -------
    x : array

    """
    if not (len(frequency)) == (len(amplitude)):
        raise TypeError("Frequency array must be same length as amplitude array")
    x = np.zeros((size))
    for j in range(len(frequency)):
        for i in range(size):
            x[i] = x[i] + amplitude[j] * np.cos(frequency[j] * i - phase)  
    return x

# def sinusoidal(size, frequency, amplitude, phase):
#     x = np.zeros((size))
#     y = np.zeros((size))
    
#     frequency_x = frequency[0]
#     frequency_y = frequency[1]
    
#     amplitude_x = amplitude[0]
#     amplitude_y = amplitude[1]
    
#     phase_x = phase[0]
#     phase_y = phase[1]
#     for i in range(0, size):
#         x[i] = amplitude_x * np.cos(frequency_x * i - phase_x) #next x-value, using the new phase
#         # x[i] = x[i] + noise_test(x[i])
#         y[i] = amplitude_y * np.sin(frequency_y * i - phase_y) #and new y-value of coord.
#         # y[i] = y[i] + noise_test(y[i])
#     return x, y


def prep_func(image, CCD, sub_pixel, wl_ran):
    spec = read_out(bin_sum(image*CCD, sub_pixel)+noise1d(bin_sum(image, sub_pixel))) 
    spec = int_r(spec, wl_ran) #interpolate to higher resolution
    return spec

def transmission_spec_func(spectrum1, spectrum2, wl_ran, disper, slitpos, img_size):
    """
    Function used to process two stellar spectrum, so it is possible to analyze the transmission spectrum of an exoplanetary
    atmosphere. 
    First, the two spectra will be summed, and read-out (2D image is collapsed column-wise). Then, the spectra are shifted
    so they have the largest correlation (best alignment). Afterwards, a linear regression is made to find the wavelength/pixel
    relation. and a moving mean filter is overlaid to smooth.

    """
    import scipy.signal
    # import input_file as inp
    
    # CCD = np.load(inp.in_CCD)
    # rout = bin_sum(image, sub_pixel)
    # r1 = read_out(rout) 
    
    # rin = bin_sum(image2, sub_pixel)
    # r2 = read_out(rin) #sum image and readout 
    # r1, r2 = int_r(r1, r2, wl_ran) #Interpolate to higher resolution

    # if noiseinp == "y":
    #     no =  noise1d(rout)
    #     ni = noise1d(rin)
    # else: 
    #     no=0
    #     ni=0   
        
          
    # if move == "y":
    autocor = scipy.signal.correlate(spectrum1, spectrum1, mode="same") #perform autocorr.
    cor = scipy.signal.correlate(spectrum1, spectrum2, mode="same") #Regular correlation
    first = np.argmax(autocor)
    second = np.argmax(cor)
    delta = first-second #amount of sub-pixels to move r1, for the two spectra to overlap
        
        # if noiseinp == "y":

    # r1, r2 = int_r(spectrum1, spectrum2, wl_ran)
    spectrum1 = np.roll(spectrum1, delta) #Move r1
        
    del first, second, autocor, cor#, rout, rin
    # if not move == "y":
    #     delta = 0

    pos = (disper[0]+slitpos[0])*100.0 #Position of each wavelength on the detector
    from scipy.stats import linregress
    a, b, r, p, s = linregress(pos, np.arange(wl_ran[0], wl_ran[1])) #Linear regression to find the lambda/pixel correlation
    wavelength = a*np.arange(img_size[1]*100.0)+(b)
    del a, b, r, p, s, 
    
    wave = wavelength[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1] #remove outlying entries, where the spectrum is not present (belo 300 nm, and above 1000) 
    spectrum1 = spectrum1[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1] 
    spectrum2 = spectrum2[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1]
    
    from astropy.convolution import convolve, Gaussian1DKernel
    spectrum1 = convolve(spectrum1,kernel = Gaussian1DKernel(4246.6)) #Moving Mean filter by convolution. Kernel is Gaussian, input is sigma

    spectrum2 = convolve(spectrum2,kernel = Gaussian1DKernel(4246.6)) #https://docs.astropy.org/en/stable/convolution/
    return spectrum1, spectrum2, wave, delta


def photon_convert(wavelength_array, flux_array, stellar_radius, distance):
    """
    Function to convert stellar flux from SPECTRUM into photon counts per second per cm^2. 

    Parameters
    ----------
    wavelength_array : array
        Array with each entry being the wavelength, in cm
    flux_array : array
        SPECTRUM fluxes. Has to be in erg/s/cm^2/Å
    stellar_radius : float
        Stellar radius of the target star. 
    distance : float
        Distance to target star, in the same unit as stellar_radius

    Returns
    -------
    spec : array
        Photon counts per second per cm^2 in the specified wavelength. [No. of photons/s/cm^2]. 
    """
    import astropy.constants as co
    
    flux = np.zeros((flux_array.shape[0]))
    for i in range(flux_array.shape[0]):
        flux[i] = np.pi*flux_array[i]*(stellar_radius/distance)**2 #The pi is a geometric factor. See 1979ApJS...40....1K
    spec = wavelength_array*flux/(co.h.cgs.value * co.c.cgs.value)
    return spec


#Trash or old functions be here:
'''
# def transmission_spec_func(image, image2, sub_pixel, wl_ran, disper, slitpos, img_size, move="y", noiseinp="n"):
#     """
#     Function used to process two stellar spectrum, so it is possible to analyze the transmission spectrum of an exoplanetary
#     atmosphere. 
#     First, the two spectra will be summed, and read-out (2D image is collapsed column-wise). Then, the spectra are shifted
#     so they have the largest correlation (best alignment). Afterwards, a linear regression is made to find the wavelength/pixel
#     relation. and a moving mean filter is overlaid to smooth.

#     """
#     import scipy.signal
#     import input_file as inp
#     CCD = np.load(inp.in_CCD)
#     rout = bin_sum(image, sub_pixel)
#     r1 = read_out(rout) 
    
#     rin = bin_sum(image2, sub_pixel)
#     r2 = read_out(rin) #sum image and readout 
#     r1, r2 = int_r(r1, r2, wl_ran) #Interpolate to higher resolution

#     if noiseinp == "y":
#         no =  noise1d(rout)
#         ni = noise1d(rin)
#     else: 
#         no=0
#         ni=0   
        
          
#     if move == "y":
#         autocor = scipy.signal.correlate(r1, r1, mode="same") #perform autocorr.
#         cor = scipy.signal.correlate(r1, r2, mode="same") #Regular correlation
#         first = np.argmax(autocor)
#         second = np.argmax(cor)
#         delta = first-second #amount of sub-pixels to move r1, for the two spectra to overlap
        
#         if noiseinp == "y":
#             rout = read_out(bin_sum(image*CCD, sub_pixel)+no) 
#             rin = read_out(bin_sum(image2*CCD, sub_pixel)+ni)
#             r1, r2 = int_r(rout, rin, wl_ran)
        
#         r1 = np.roll(r1, delta) #Move r1
        
#         del first, second, autocor, cor, rout, rin
#     if not move == "y":
#         delta = 0

#     pos = (disper[0]+slitpos[0])*100.0 #Position of each wavelength on the detector
#     from scipy.stats import linregress
#     a, b, r, p, s = linregress(pos, np.arange(wl_ran[0], wl_ran[1])) #Linear regression to find the lambda/pixel correlation
#     wavelength = a*np.arange(img_size[1]*100.0)+(b)
#     del a, b, r, p, s, 
    
#     wave = wavelength[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1] #remove outlying entries, where the spectrum is not present (belo 300 nm, and above 1000) 
#     r1 = r1[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1] 
#     r2 = r2[np.max(np.where(wavelength<wl_ran[0]))+1:np.min(np.where(wavelength>wl_ran[1]))-1]
    
#     from astropy.convolution import convolve, Gaussian1DKernel
#     r1 = convolve(r1,kernel = Gaussian1DKernel(4246.6)) #Moving Mean filter by convolution. Kernel is Gaussian, input is sigma

#     r2 = convolve(r2,kernel = Gaussian1DKernel(4246.6)) #https://docs.astropy.org/en/stable/convolution/
#     return r1, r2, wave, delta




def disperser(wl_endpoints, jit_img, psf_img, pos, image_size, dispersion, eff, 
              mask_img, steps=1, secondary_source='n', plot='n'):
    import sys
    from scipy import signal
    x_pos=pos[0] 
    y_pos=pos[1] #load in position of "zeroth order"
    im_disp = np.zeros((image_size[0],image_size[1])) # empty image
    im_disp_lambda = np.zeros((image_size[0],image_size[1])) 
    x_dispersion = dispersion[0] #load in dispersions
    y_dispersion = dispersion[1]
    numColors = int( (wl_endpoints[1]-wl_endpoints[0])) #total number of colours to iterate
    print("Number of colors to iterate: " + str(numColors))
    print(' ')
    if plot=='y': #this part is not useful atm
        import matplotlib.pyplot as plt
        plt.figure()
        from matplotlib.colors import LinearSegmentedColormap
        N = 256 #8-bit value, to fix colours
        colspec = plt.cm.get_cmap('Spectral') #Fetches colourmap to use later
        vals = np.ones((N,4)) #Setup for colormap

    for i in range(0, numColors, steps):
    # for i in range(0,101, steps):
        im = np.zeros((image_size[0],image_size[1])) #create temp. image
        if secondary_source == 'y': #To account for the secondary light source perhaps not being fully within the psf
            # fold = folding(psf_img[:,:,i], jit_img)
            fold = signal.convolve2d(psf_img[:,:,i], jit_img, mode='same', boundary='fill') #fold psf and jitter
            fold = fold[0:jit_img.shape[1], 0:jit_img.shape[0]] #cut down to regular shape
        else:
            fold = signal.convolve2d(psf_img[:,:,i], jit_img, mode='same', boundary='fill') #fold as usual, if no sec. sources
        # fold=fold/np.sum(fold)
        
        foo = int(psf_img.shape[0]/2)
        # im[0+x_pos-foo:len(jitter)+x_pos-foo, 0+y_pos-foo:len(jitter)+y_pos-foo] = im[0+x_pos-foo:len(jitter)+x_pos-foo, 0+y_pos-foo:len(jitter)+y_pos-foo] + fold*magni
        im[0+y_pos-foo:len(fold)+y_pos-foo, 0+x_pos-foo:len(fold)+x_pos-foo] = fold #im[0+y_pos-foo:len(fold)+y_pos-foo, 0+x_pos-foo:len(fold)+x_pos-foo] + fold#*magni
        immask = im*mask_img #mask is "overlaid" by multiplying
        
        roll_x = np.roll(immask,  int(np.modf(x_dispersion[i])[1]), axis=1) #move/disperse the light 
        roll_y = np.roll(roll_x,  int(np.modf(y_dispersion[i])[1]), axis=0) #also in the y-direction

        dx = abs(np.modf(x_dispersion[i])[0]) #residual amount (decimal amounts are shifted to the next sub-pixel)
        dy = abs(np.modf(y_dispersion[i])[0])
        
        foob = roll_y*(eff[i]*(1-dx)*(1-dy)) #multiply by efficiency 
        im_disp = im_disp + foob  # Add the rolled image to the final, and multiply by the "effectivity"
        
        
        roll_dx = np.roll(roll_y, 1, axis=1) # Roll the residual to the next subpixel
        eff_dx = eff[i] * dx * (1-dy) # effectivity of the x-residual
        
        roll_dy = np.roll(roll_y, 1, axis=0) # Roll the residual to the next subpixel, y-wise
        eff_dy = eff[i] * dy * (1-dx) # y-residual eff.
        
        roll_dxy = np.roll(roll_dx, 1, axis=0) # roll the image one step in both x- and y-wise.
        eff_dxy = eff[i]* dx * dy   #and eff.
        
        baar = roll_dx*eff_dx + roll_dy*eff_dy + roll_dxy*eff_dxy
        
        im_disp = im_disp + baar #add all residuals and multiply by their respective effectivities.
        
        im_disp_lambda = im_disp_lambda+((foob+baar)*(i+wl_endpoints[0])) #fill in im_disp, and multiply by wavelength i
        # im_disp_lambda = im_disp_lambda+(i+wl_endpoints[0]) #fill in im_disp, and multiply by wavelength i
        sys.stdout.write('/'); sys.stdout.flush(); #"Progress bar", just for visuals    
        
        ##### Plotting #####
        if plot == 'y':
            vals[:, 0] = np.linspace(0, colspec(1-i/750)[0], N) #Making new colourmap values
            vals[:, 1] = np.linspace(0, colspec(1-i/750)[1], N) #the /750 is to normalize the colormap, so values fall between 0 and 1
            vals[:, 2] = np.linspace(0, colspec(1-i/750)[2], N)
            vals[:, 3] = np.linspace(0, 1, N) #alpha, for making the cmap transparent
            newcmp = LinearSegmentedColormap.from_list(name='Spectral', colors=vals) #Creates new cmp, based on vals    
            plt.imshow(roll_y, cmap=newcmp) # Show array   
    
    if plot=='y':
        plt.title('Color dispersion of sample spectrum', size=18)
        plt.xlabel('Sub-pixel', size=13)
        plt.ylabel('Sub-pixel', size=13)
    return im_disp, im_disp_lambda

'''


"""    
def jitter(steps=1000, dt=10, time_delay=5, gain=0.2, amplitude_act=0.1, amplitude_sens=0.1):
    Jitter generator
    
    Parameters
    ----------
    steps : int
        desired number of entries in position vector
    dt : int
        "batch size", how many indices to run through before correction
    time_delay : int
        index that is subtraced from the dt'th index, used to correct for. 
    gain : float
        gain of correction. RoT - should be around 1/time_delay
    amplitude_* : float
        size of noise added under correction
    
    Returns
    -------
    x, y : array_like
        Vectors with x and y positions, of the size specified (+1)
    
    
    x = np.zeros(steps+1) #Allocates  vectors for x and y position
    y = np.zeros(steps+1)
    k = 0
    for j in range(int(steps/dt)):
        jitt = np.random.randn(1,2) #Generate random noise to add to position
        for i in range(1,dt): #Jitter will be added to position, cumulatively
            x[k+i] = x[k+i-1]+jitt[0,0] #Takes previous position, adds jitter
            y[k+i] = y[k+i-1]+jitt[0,1]
            jitt = np.random.randn(1,2)*0.05 #Generate new jitter, to add to next postition
        x_correction = gain*(-x[k+i-time_delay])+amplitude_act*np.random.randn()+amplitude_sens*np.random.randn() #Generates the correction term, 
        # but for the index "time_delay" ago - so for time_delay= 5, x[k+9-5] = x[k+4].  
#        print(x_correction)
        y_correction = gain*(-y[k+i-time_delay])+amplitude_act*np.random.randn()+amplitude_sens*np.random.randn()
            
        x[k+i+1] = x[k+i] + x_correction #correction term is added to the last entry in the small batch of "steps"
        y[k+i+1] = y[k+i] + y_correction
        k=k+dt #K is updated, and the whole thing runs again, this time for index +dt. 
    x = x[0:steps]
    y = y[0:steps] #Cut off the last step, as it is just filler
    return x, y

def PSF(sigma_x, sigma_y, res=100, x_size=30, y_size=30, x0=0, y0=0):
    """ """This function creates a mesh of the PSF as a 2D Gaussian
    Parameters
    ----------
    sigma_x, sigma_y : float
        std. The "spread" of the function
    res : int
        number of steps in grid. "Resolution" of PSF
    x_size, y_size : int
        Size of created grid
    x0, y0 : float
        x and y position on the detector
    """ """
    import numpy as np
    x = np.linspace(-x_size, x_size, res)
    y = np.linspace(-y_size, y_size, res)
    x, y = np.meshgrid(x, y) #define meshgrid
    z = (1/(2*np.pi*sigma_x*sigma_y) * np.exp(-((x)**2/(2*sigma_x**2) 
         + (y)**2/(2*sigma_y**2)))) #2 Gaussian functions
    z = z/np.sum(z)
    return x, y, z


def psf_maker(res=101, wl_start=350, wl_stop=1100, si=[30, 30]): 
    """ """
    Creates 10 psf's. Mainly for quickly testing the psf-interpolator.

    Parameters
    ----------
    res : int, optional
        Size of the psf. The default is 101.
    wl_start : int, optional
        Start wavelength. The default is 350.
    wl_stop : int, optional
        Stop wavelength. The default is 1100.
    si : TYPE, optional
        DESCRIPTION. The default is [30, 30].

    Returns
    -------
    psf_temp : array
        PSF to interpolate.
    ran : array
        List of corresponding wavelengths. Each entry is the wavelength of the PSF slice in psf_temp
    """ """
    #Load packages
    import numpy as np

    ran = np.linspace(wl_start, wl_stop, 10) #set for loop index
    foo = np.zeros((res, res)) #create empty grid to store in the psf
    for i in ran:
        sigma_x = np.log(i)+0.5*i/100 # Used in the 2D Gaussian
        sigma_y = np.log(i)+0.5*i/100
        x, y, z = PSF(sigma_x, sigma_y, res=res, x_size=si[0], y_size=si[1]) #2D Gaussian
        foo = np.dstack((foo,z)) #appends the new psf to the 3d matrix
    psf_temp = foo[:,:,1:ran.shape[0]+1] #cuts off the initial layer of 0's
    del foo    
    print(' ')
    print('psf done')
    return  psf_temp, ran
    
    
def disp_func(wl_start, wl_stop, delta_lambda):
    return 1.15*np.arange(wl_start, wl_stop, delta_lambda)-800


def disperser(image, psf, magni, mask, eff, dispers=(350, 1100, 1), CCDsize=(1000,1000), cols=(0, 750), stepsize=5, exposure=100, plot='n', save='n'):
    import sys
    if plot == 'y':
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
        N = 256 #8-bit value, to fix colours
        colspec = plt.cm.get_cmap('Spectral') #Fetches colourmap to use later
        vals = np.ones((N,4))
        plt.figure()
    
    im_disp = np.zeros(CCDsize)
    x_j, y_j = jitter(steps=100, gain=0.2, amplitude_act=3, amplitude_sens=3)
    ux = int(np.floor(psf.shape[0]/2))
    ox = int(np.floor(psf.shape[0]/2)+1)
    uy = int(np.floor(psf.shape[1]/2))
    oy = int(np.floor(psf.shape[1]/2)+1)
    
    wl_start = dispers[0]
    wl_stop = dispers[1]
    delta_lambda = dispers[2]
    dispersion=disp_func(wl_start, wl_stop, delta_lambda)
    print(' ')
    print('Number of colours complete:')
    print('         10        20        30        40        50        60        70        80        90')
    for k in range(cols[0],cols[1], stepsize):
        x_pos = 499 #generates star position
        y_pos = 499 #has to be here, otherwise jitter won't be applied properly
        x_0 = x_pos
        y_0 = y_pos
        for i in range(exposure):
            image[x_pos-ux:x_pos+ox, y_pos-uy:y_pos+oy] = image[x_pos-ux:x_pos+ox, y_pos-uy:y_pos+oy]+psf[:,:,k]*magni #adds psf values to selected area of image array
            x_pos = x_0+x_j[i]
            y_pos = y_0+y_j[i] #updates coordinates based on jitter
            x_pos = int(np.around(x_pos))
            y_pos = int(np.around(y_pos)) # rounds off the coordinates, as matrix can only take int as index
        image_masked = image[:,:]*mask #Overlay slit mask
        roll = np.roll(image_masked, int(dispersion[k]), axis=1)*eff[k]
        im_disp = im_disp + roll + np.random.standard_normal((1000, 1000))*0.001 #Disperses the colours, using np.roll
        
        sys.stdout.write('/'); sys.stdout.flush(); #"Progress bar", just for visuals    
        
        ##### Plotting #####
        if plot == 'y':
            vals[:, 0] = np.linspace(0, colspec(1-k/750)[0], N) #Making new colourmap values
            vals[:, 1] = np.linspace(0, colspec(1-k/750)[1], N) #the /750 is to normalize the colormap, so values fall between 0 and 1
            vals[:, 2] = np.linspace(0, colspec(1-k/750)[2], N)
            vals[:, 3] = np.linspace(0, 1, N) #alpha, for making the cmap transparent
            newcmp = LinearSegmentedColormap.from_list(name='Spectral', colors=vals) #Creates new cmp, based on vals    
            plt.imshow(roll, cmap=newcmp) # Show array
        
    if plot == 'y':
        plt.title('Color dispersion of sample spectrum', size=18)
        plt.xlabel('Sub-pixel', size=13)
        plt.ylabel('Sub-pixel', size=13)
    if save == 'y':
        plt.savefig('disp.png', dpi=400)    
    return im_disp
    
def disperser2_copy(jit_img, psf_img, pos, image_size, dispersion, eff, magni, mask_img, steps=1, plot='n'):
    '''
    Parameters
    ----------
    jitter : array of float64
        Jitter "image". 
    psf : _hl.dataset.Dataset
        Point Spread Function image. Must be a 3D array with depth equal to the number of colors
    pos : list
        Position of the star. Two values in a list.
    image_size : tuble or list
        Two integers that determine the size of the image. Must have the same dimensions as the mask and jitter
    dispersion : tuble
        Requires two entries, one for dispersion in the x- and one for the y-direction. Must have same length as number of colors
    eff : array of float64
        Spectral effeciency/throughput. Must be same lenght as number of colors
    magni : float
        Magnitude of the star.
    mask : array of float64
        Slit mask. Must have same dimensions as image.
    steps : int, optional
        Size of color "steps" to include in the disperser. The default is 1 - so all colors are included.
    plot : string, optional
        Toggles plotting of the color-dispersion, mainly for visuals. The default is 'n'.

    Returns
    -------
    im_disp : array of float64
        Dispersed image, 2D array.

    '''
    import sys
    x_pos=pos[0]
    y_pos=pos[1]
    im_disp = np.zeros((image_size[0],image_size[1]))
    x_dispersion = dispersion[0]
    y_dispersion = dispersion[1]
    
    if plot=='y':
        import matplotlib.pyplot as plt
        plt.figure()
        from matplotlib.colors import LinearSegmentedColormap
        N = 256 #8-bit value, to fix colours
        colspec = plt.cm.get_cmap('Spectral') #Fetches colourmap to use later
        vals = np.ones((N,4)) #Setup for colormap

    for i in range(0,750, steps):
    # for i in range(0,101, steps):
        im = np.zeros((image_size[0],image_size[1]))
        fold = folding(psf_img[:,:,i], jit_img)
        fold=fold/np.sum(fold)
        
        foo = int(psf_img.shape[0]/2)
        # im[0+x_pos-foo:len(jitter)+x_pos-foo, 0+y_pos-foo:len(jitter)+y_pos-foo] = im[0+x_pos-foo:len(jitter)+x_pos-foo, 0+y_pos-foo:len(jitter)+y_pos-foo] + fold*magni
        im[0+y_pos-foo:len(fold)+y_pos-foo, 0+x_pos-foo:len(fold)+x_pos-foo] = fold #im[0+y_pos-foo:len(fold)+y_pos-foo, 0+x_pos-foo:len(fold)+x_pos-foo] + fold#*magni
        immask = im*mask_img
        
        roll_x = np.roll(immask, int(np.modf(  x_dispersion[i])[1]), axis=1)
        roll_y = np.roll(roll_x,  int(np.modf(y_dispersion[i])[1]), axis=0)

        dx = abs(np.modf(x_dispersion[i])[0])
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
        
        ##### Plotting #####
        if plot == 'y':
            vals[:, 0] = np.linspace(0, colspec(1-i/750)[0], N) #Making new colourmap values
            vals[:, 1] = np.linspace(0, colspec(1-i/750)[1], N) #the /750 is to normalize the colormap, so values fall between 0 and 1
            vals[:, 2] = np.linspace(0, colspec(1-i/750)[2], N)
            vals[:, 3] = np.linspace(0, 1, N) #alpha, for making the cmap transparent
            newcmp = LinearSegmentedColormap.from_list(name='Spectral', colors=vals) #Creates new cmp, based on vals    
            plt.imshow(roll_y, cmap=newcmp) # Show array   
    
    if plot=='y':
        plt.title('Color dispersion of sample spectrum', size=18)
        plt.xlabel('Sub-pixel', size=13)
        plt.ylabel('Sub-pixel', size=13)
    return im_disp
"""