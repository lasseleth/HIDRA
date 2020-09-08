import numpy as np
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

def BMB_interp(x, y, wl_start=300, wl_stop=1200, delta_lambda=0.01, kind='cubic', lowlim=400, uplim=1100):
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
    wl_start : float
        shortest wavelength, in nm.
    wl_stop : float
        longest wavelength, in nm.
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
    
    Examples
    --------
    >>> x = np.array([100, 217, 1000])
    >>> y = np.array([100, 217, 1000])
    >>>  BMB_interp(x, y, 100, 1000, 0.01, kind='linear', lowlim=120, uplim=950)
    array([[ 100.        ,    0.        ],
           [ 100.01000011,    0.        ],
           ...,
           [ 999.98999989,    0.        ],
           [1000.        ,    0.        ]])
    """
    from scipy.interpolate import interp1d #Load neccessary package
    import numpy as np
    f = interp1d(x, y, kind=kind, fill_value="extrapolate") #interpolates, and extrapolates if the given table does not cover the wavelength range
    xnew = np.linspace(wl_start, wl_stop, num=int((wl_stop-wl_start)/delta_lambda), endpoint=True) #Generates new x-values
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
    
    noise = np.random.standard_normal((x_size, y_size))*var #Noise matrix 
    CCD = np.ones((x_size,y_size)) #"Clean" matrix
    CCD = CCD-noise #Subtracts noise from "clean" CCD matrix
    CCD2=CCD
    
    #Correlate the subpixels
    N = 3 # number of times to correlate
    for t in np.arange(0,N):
        for i in np.arange(0,x_size):
            for j in np.arange(0,y_size):
                pixel = CCD[i:i+S, j:j+S-1] #cuts out a pixel
                CCD2[i, j] = np.sum(np.sum(pixel)/np.size(pixel)) #correlates surrounding subpixels
        sys.stdout.write('.'); sys.stdout.flush(); #"Progress bar", just for visuals   
    #Introduces grid, to mimic the actual pixels - seperate the subpixels by a grid with a slight loss defined by grid_loss variable. 
    grid = np.ones((CCD.shape[0], CCD.shape[1])) #Set up grid
    grid[0::gridsize,:]=grid_loss #Sets gridloss for every 'gridsize' row (10)
    grid[:,0::gridsize]=grid_loss #sets gridloss for every 'gridsize' coloumn (10)
    #to see a visualization of this, use the variable explorer - type: %varexp --imshow grid
    noise2 = np.random.standard_normal((x_size, y_size))*var2 
    CCD2 = CCD2+noise2+1
    # CCD2 = CCD2/np.mean(CCD2)
    CCD = CCD2*grid #overlays the grid on the CCD
    CCD = CCD/np.mean(CCD)
    return CCD

def PSF(sigma_x, sigma_y, res=100, x_size=30, y_size=30, x0=0, y0=0):
    """This function creates a mesh of the PSF as a 2D Gaussian
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
    """
    import numpy as np
    x = np.linspace(-x_size, x_size, res)
    y = np.linspace(-y_size, y_size, res)
    x, y = np.meshgrid(x, y) #define meshgrid
    z = (1/(2*np.pi*sigma_x*sigma_y) * np.exp(-((x)**2/(2*sigma_x**2) 
         + (y)**2/(2*sigma_y**2)))) #2 Gaussian functions
    z = z/np.sum(z)
    return x, y, z


def psf_maker(res=101, wl_start=350, wl_stop=1100, si=[30, 30]): 
    """
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
    """
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

def PSF2(file_name, wl_endpoints=(350, 1100), res=101, size=(100, 100), step=1):
    """
    Creates a new file containing the full-color PSF, without interpolating as with psf_maker

    Parameters
    ----------
    file_name : str
        Desired name of the file.
    wl_endpoints : tuple, optional
        Two values that mark the first and last colors. The default is (350, 1100).
    res : int, optional
        Resolution of the PSF. The default is 100.
    size : tuple, optional
        Size of the meshgrid used in the 2D Gaussian. Just a tweakable parameter. The default is (100, 100).
    step : int, optional
        Number of color steps. The default is 1.

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
    
    ran = range(wl_endpoints[0], wl_endpoints[1], delta_lambda) #set for-loop range
    res = input_psf_images.shape[0] # Width of the input psf, so the created psf will have the same size

    psf = np.zeros((input_psf_images.shape[0], input_psf_images.shape[1], wl_endpoints[1]-wl_endpoints[0])) #Creates empty array for the new psf
    for i in range(res):
        for j in range(res):
            f_test = interp1d(input_psf_wl, input_psf_images[i,j,:], kind='quadratic') #sets up interpolation function
            psf[i,j,:] = f_test(ran) # interpolates at the wavelengths specified in the range
        sys.stdout.write('.'); sys.stdout.flush(); #"Progress bar", just for visuals
    print(' ')
    print('Interpolation done')
    print(' ')    
    return psf


def jitter(steps=1000, dt=10, time_delay=5, gain=0.2, amplitude_act=0.1, amplitude_sens=0.1):
    """Jitter generator
    
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
    """
    
    x = np.zeros(steps+1) #Allocates  vectors for x and y position
    y = np.zeros(steps+1)
    k = 0
    for j in range(int(steps/dt)):
        jitt = np.random.randn(1,2) #Generate random noise to add to position
        for i in range(1,dt): #Jitter will be added to position, cumulatively
            x[k+i] = x[k+i-1]+jitt[0,0] #Takes previous position, adds jitter
            y[k+i] = y[k+i-1]+jitt[0,1]
            jitt = np.random.randn(1,2) #Generate new jitter, to add to next postition
        x_correction = gain*(-x[k+i-time_delay]+amplitude_act*np.random.randn())+amplitude_sens*np.random.randn() #Generates the correction term, 
        # but for the index "time_delay" ago - so for time_delay= 5, x[k+9-5] = x[k+4].  
#        print(x_correction)
        y_correction = gain*(-y[k+i-time_delay]+amplitude_act*np.random.randn())+amplitude_sens*np.random.randn()
            
        x[k+i+1] = x[k+i] + x_correction #correction term is added to the last entry in the small batch of "steps"
        y[k+i+1] = y[k+i] + y_correction
        k=k+dt #K is updated, and the whole thing runs again, this time for index +dt. 
    x = x[0:steps]
    y = y[0:steps] #Cut off the last step, as it is just filler
    return x, y

def slit(slit_size=[10,100], pos=[499, 499], image_size=[1000,1000]):
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

def disp_func(wl_start, wl_stop, delta_lambda):
    return 1.15*np.arange(wl_start, wl_stop, delta_lambda)-800

"""
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
"""

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
    for i in range(len(x)):
        jitter[x[i].astype(int)+int(np.floor(psf_size[0]/2)), y[i].astype(int)+int(np.floor(psf_size[1]/2))]= jitter[x[i].astype(int)+int(np.floor(psf_size[0]/2)), y[i].astype(int)+int(np.floor(psf_size[1]/2))]+1 # Create jitter "image". +1 to every point where the jitter "hits"
    return jitter

def folding(psf_image, jitter_image):
    from scipy import signal
    folded=signal.convolve2d(psf_image, jitter_image, mode='same', boundary='fill') #convolves the psf slice and the jitter image
    return folded

#The correct disperser::::
def disperser2(wl_endpoints, jit_img, psf_img, pos, image_size, dispersion, eff, magni, mask_img, steps=1, plot='n'):
    import sys
    x_pos=pos[0] 
    y_pos=pos[1]
    im_disp = np.zeros((image_size[0],image_size[1]))
    im_disp_lambda = np.zeros((image_size[0],image_size[1]))
    x_dispersion = dispersion[0]
    y_dispersion = dispersion[1]
    numColors = int( (wl_endpoints[1]-wl_endpoints[0]))
    print(numColors)
    print(' ')
    if plot=='y':
        import matplotlib.pyplot as plt
        plt.figure()
        from matplotlib.colors import LinearSegmentedColormap
        N = 256 #8-bit value, to fix colours
        colspec = plt.cm.get_cmap('Spectral') #Fetches colourmap to use later
        vals = np.ones((N,4)) #Setup for colormap

    for i in range(0, numColors, steps):
    # for i in range(0,101, steps):
        im = np.zeros((image_size[0],image_size[1]))
        fold = folding(psf_img[:,:,i], jit_img)
        # fold=fold/np.sum(fold)
        
        foo = int(psf_img.shape[0]/2)
        # im[0+x_pos-foo:len(jitter)+x_pos-foo, 0+y_pos-foo:len(jitter)+y_pos-foo] = im[0+x_pos-foo:len(jitter)+x_pos-foo, 0+y_pos-foo:len(jitter)+y_pos-foo] + fold*magni
        im[0+y_pos-foo:len(fold)+y_pos-foo, 0+x_pos-foo:len(fold)+x_pos-foo] = fold #im[0+y_pos-foo:len(fold)+y_pos-foo, 0+x_pos-foo:len(fold)+x_pos-foo] + fold#*magni
        immask = im*mask_img
        
        roll_x = np.roll(immask,  int(np.modf(x_dispersion[i])[1]), axis=1) #move/disperse the light 
        roll_y = np.roll(roll_x,  int(np.modf(y_dispersion[i])[1]), axis=0) #also in the y-direction

        dx = abs(np.modf(x_dispersion[i])[0]) #residual amount
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
    import numpy as np
    counts = np.array(())
    for i in range(dispersed.shape[0]):
        counts = np.append(counts, np.sum(dispersed[:,i]))
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

    for x in range(0, inp.shape[0], bin_size): #Range for 1st
        j = range(0+x, bin_size+x) #Bin range. ex. 20-30 if bin_size is 10
        for i in range(0, inp.shape[0]): # over all columns 
            temp[i, int(j[0]/bin_size)]= sum(inp[i,j]) #sum, and add to temp
    
    for x in range(0, inp.shape[1], bin_size): #2nd step, repeat 1st step, but for rows
        i = range(0+x, bin_size+x) #row bin-range. 
        for j in range(0, summed.shape[0]): 
            summed[int(i[0]/bin_size), j]= sum(temp[i,j]) #sum and add to result-matrix
            
    return summed
