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

