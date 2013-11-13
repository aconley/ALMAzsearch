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
    n3, sn3 = twoline.sn(tuning3, z, lir=lir, applyciratio=True)

Now n3 shows the number of lines detected at each redshift (a maximum of 2),
and sn3 gives the S/N ratios for each line.

It is also possible to specify a target 500um flux density instead
of a L_IR using the snu argument to `twoline.sn`.  This then uses
a dust model based on typical DSFG parameters for red sources
(see Dowell et al. 2013), although these can be adjusted.

### Dependencies
This depends on a number of python packages:

* [numpy](http://www.numpy.org/) version 1.7 or later.
* [scipy](http://www.scipy.org/)
* [astropy](http://www.astropy.org/) version 0.3 or later.
* [cython](http://cython.org/)

### References
* More information about ALMA can be found at the
  [ALMA Project Pages](http://www.almascience.nrao.edu)
