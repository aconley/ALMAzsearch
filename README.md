ALMAzsearch
===========

Code to plan CO redshift searches in spectral scan mode
with ALMA in Cycle 2 or later.

###Installation
The usual

        python setup.py install

###Usage
First, you set up tunings, which are setups for ALMA spectral scans.

    import astropy.units as u
    import ALMAzsearch
    import numpy as np
    tuning3 = ALMAzsearch.scan_tuning(u.Quantity(84, u.GHz),
                                      u.Quantity(5, u.min))	

Then you chose your line templates, and set up and execute 2-line search
     
    line_templates = ["Eyelash", "FLS3", "HLSW01"] # Names of line templates
    twoline = ALMAzsearch.twoline(line_templates, ciratio=1.5)
    z = np.linspace(3.5, 7.0, 500) # Redshift range to search
    lir = u.Quantity(3e13, u.solLum)
    n3, sn3 = twoline.sn(tuning3, z, norm=lir, applyciratio=True)

Now n3 shows the number of lines detected at each redshift (at any S/N)
and sn3 gives the S/N ratios for each line as a nz by 2 array.
The line list includes 12CO and [CI] lines only.

It is also possible to specify a target 500um (or other wavelength)
flux density instead of a L_IR using the norm argument to `twoline.sn`.  
This then uses a dust model based on typical DSFG parameters for red sources
(see Dowell et al. 2013), although these can be adjusted.  For example,
to scale to a 15 mJy SCUBA2 (850um) flux density:

    obs_wave = u.Quantity(850, u.um)
    snu = u.Quantity(15., u.mJy)
    n3, sn3 = twoline.sn(tuning3, z, norm=snu, obs_wave=obs_wave,
                         applyciratio=True)

    
### Dependencies
This depends on a number of python packages:

* [numpy](http://www.numpy.org/) version 1.7 or later.
* [scipy](http://www.scipy.org/)
* [astropy](http://www.astropy.org/) version 0.3 or later.
* [cython](http://cython.org/)

### References
* More information about ALMA can be found at the
  [ALMA Project Pages](http://www.almascience.nrao.edu)
