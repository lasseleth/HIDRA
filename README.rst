=====
HIDRA
===== 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Hyperspectral Instrument Data Resemblance Algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This python 3.6+ package is capable of simulation spectroscopic observations of stellar spectra, from photons arriving at the telescope, to the photo-electrons being detected by a CCD. 

The scripts were created as a part of my Master's thesis, which was titled "Simulating multi-wavelength observations from low-resolution spectrographs".

Author: Lasse L. S. Berthelsen

- `More documentation <https://hidra.readthedocs.io/en/latest/>`_

Installation
~~~~~~~~~~~~

The best way to install ``HIDRA`` is

.. code-block:: bash

   git clone https://github.com/lasseleth/HIDRA.git
   cd Sim

Simulating a stellar observation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Given a 1D spectrum of a target star, this package produces a 2D 
image of how this spectrum would appear on a CCD. 
To make a simple simulation, my code might look like:

.. code-block:: python

   # Imports
   import HIDRA
   import input_file as inp
   
   #Setup the code using user-defined inputs. These can be adjusted in the input file
   spec_eff, spec_eff2, jitter, x_j, y_j, psf, img_size, sub_pixel, pl_arc_pix, disper, mask, slitpos, background = HIDRA.setup(inp)
   
   image, image_wl=HIDRA.disperser(wl_endpoints=wl_ran, jit_img=jitter, psf_img=psf, pos=slitpos, image_size=img_size, 
                                        dispersion=disper, eff=spec_eff, mask_img=mask, steps=1, plot='n')
