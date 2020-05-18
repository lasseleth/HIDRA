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
    This function creates a CCD composed of subpixels, with a separating grid between all full pixels. The grid will have some loss.
    
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
    
    CCD = CCD2*grid #overlays the grid on the CCD
    
    return CCD

def PSF(sigma_x=1, sigma_y=1, res=100, x_size=30, y_size=30, x0=0, y0=0):
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
    return x, y, z


def psf_maker(file_path="/home/lasse/Documents/uni/Speciale/Sim/psf_array.hdf5", res=101, wl_start=350, wl_stop=1100, delta_lambda=1):    
    #Load packages
    import h5py #package for handling the upcoming HUGE arrays, that the RAM won't be able to handle
    import numpy as np
    import sys
    from scipy.interpolate import interp1d
    import os
    
    if os.path.exists(file_path) == True: #If file already exists, it will be deleted
        os.remove(file_path)
    
    #and then remade
    psf_file = h5py.File(file_path, "a")
    psf_file.create_dataset('psf', ( res, res, int((wl_stop-wl_start)/delta_lambda) ), dtype='f') #creates datasets in the file. 
    #This is for the PSF generated a bit later
    
    psf   = psf_file['psf']

    ran = np.linspace(wl_start, wl_stop, 10) #set for loop index
    foo = np.zeros((res, res)) #create empty grid to store in the psf
    for i in ran:
        x, y, z = PSF(sigma_x=np.log(i+2), sigma_y=np.log(i+2), res=res) #2D Gaussian
        foo = np.dstack((foo,z)) #appends the new psf to the 3d matrix
    psf_temp = foo[:,:,1:ran.shape[0]+1] #cuts off the initial layer of 0's
    del foo


    for i in range(res):
        for j in range(res):
            f_test = interp1d(ran/delta_lambda, psf_temp[i,j,:], kind='quadratic') #sets up interpolation
            psf[i,j,:] = f_test(range(int(wl_start/delta_lambda), int(wl_stop/delta_lambda))) # interpolates at the wavelengths specified in the range
        sys.stdout.write('.'); sys.stdout.flush(); #"Progress bar", just for visuals
    
    print(' ')
    print('psf done')
    del psf_temp, f_test, x, y, z, ran
    print(psf_file['psf'][0,0,0])
    
    return psf_file
    