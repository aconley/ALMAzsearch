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
       tuning3 = ALMAzsearch.scan_tuning(u.Quantity(84, u.GHz),
                                  	 u.Quantity(5, u.min))	

Then you chose your line templates, and set up and execute 2-line search
     
        line_templates = ["Eyelash", "FLS3", "HLSW01"]
	twoline = ALMAzsearch.twoline(line_templates, ciratio=1.5)
	z = np.arange(3.5, 7.0, 500)
	n3, sn3 = twoline.sn(tuning3, z, lir=lir, applyciratio=True)

Now n3 shows the number of lines detected at each redshift (a maximum of 2),
and sn3 gives the S/N ratios for each line.

### Dependencies
This depends on a number of python packages:
* [numpy](http://www.numpy.org/) version 1.7 or later.
* [scipy](http://www.scipy.org/)
* [astropy](http://www.astropy.org/) version 0.3 or later.

### References
* More information about ALMA can be found at the
  [ALMA Project Pages](http://www.almascience.nrao.edu)
