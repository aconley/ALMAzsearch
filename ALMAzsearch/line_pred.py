import numpy as np
import astropy.units as u
from astropy.cosmology import WMAP9

"""Predicts line intensities at a given redshift from a template"""

from .co_lines import linefreq

__all__ = ["line_template"]

class line_template(object) :
    def __init__(self,name,ciratio=1.0) :
        # We store which lines are observed (or predicted) as
        # well as L_IR and the observed line strengths (both
        # corrected for lensing)
        if name == "Arp220" :
            self.source = "Rangwala et al. 2011"
            self.z = 0.0181
            self.dl = u.Quantity(77, u.Mpc)
            self.lir = u.Quantity(1.77e12, u.solLum)
            
            #Names of lines
            self.linename = ['12CO(1-0)','12CO(2-1)','12CO(3-2)',
                             '12CO(4-3)','12CO(5-4)','12CO(6-5)',
                             '12CO(7-6)','12CO(8-7)','12CO(9-8)',
                             '12CO(11-10)','12CO(12-11)','12CO(13-12)',
                             '[CI](1-0)','[CI](2-1)']


            #W m^-2, observed
            self.linei0 = u.Quantity([1.582e-18,  8.512e-18,  4.162e-17,
                                      6.9895e-17, 7.0314e-17, 9.3966e-17,
                                      9.3205e-17, 10.014e-17, 10.1133e-17,
                                      5.313e-17,  3.803e-17,  2.207e-17,
                                      3.015e-17, 5.73e-17], u.W / u.m**2)
            

            self.lineratio = numpy.ones_like(self.linei0)
            self.lineratio[-2:] = ciratio

        elif name == "Eyelash" :
            self.source = "Danielson et al. 2011"
            self.z = 2.3259
            self.dl = WMAP9.luminosity_distance(self.z)
            self.lir = u.Quantity(2.3e12, u.solLum)

            #Names of lines
            self.linename = ['12CO(1-0)','12CO(3-2)','12CO(4-3)',
                             '12CO(5-4)','12CO(6-5)','12CO(7-6)',
                             '12CO(8-7)','12CO(9-8)',
                             '[CI](1-0)','[CI](2-1)']

            #W m^-2, observed, corrected for lensing
            self.linei0 = u.Quantity([7.68e-23, 1.409e-21, 2.461e-21, 3.323e-21,
                                      4.58e-21, 3.13e-21, 2.5e-21, 1.15e-21,
                                      2.43e-21, 4.044e-21], u.W / u.m**2)
            self.lineratio = numpy.ones_like(self.linei0)
            self.lineratio[-2:] = ciratio

        elif name == "HLSW01" :
            self.source = "K.S. Scott et al. 2011"
            self.z = 2.9575
            self.lir = u.Quantity(1.43e13, u.solLum)
            self.dl = WMAP9.luminosity_distance(self.z)

            #2-1,4-3,6-5,8-7 are model based
            #[CI](1-0) assumes 30% of 4-3
            self.linename = ['12CO(1-0)','12CO(2-1)','12CO(3-2)',
                             '12CO(4-3)','12CO(5-4)','12CO(6-5)',
                             '12CO(7-6)','12CO(8-7)','12CO(9-8)',
                             '12CO(10-9)','[CI](1-0)']
            self.linei0 = u.Quantity([1e-22,   1e-21,    2.6e-21,
                                      0.6e-20, 1.05e-20, 1.5e-20,
                                      2.2e-20, 1.7e-20,  1e-20, 
                                      1.3e-20, 0.2e-20], u.W / u.m**2)
            self.lineratio = numpy.ones_like(self.linei0)
            self.lineratio[-1:] = ciratio

        elif name == "ID141" :
            #Observer frame SED params; note none of the below are
            # corrected for lensing.  The CO (6-5) and (3-2) are from
            # their model.  Since the lensing magnification is not known,
            # a value of mu = 5 has been assumed
            mu = 5.0
            self.source = "Cox et al. 2011"
            self.z = 4.243
            self.lir = u.Quantity(8.5e13 / mu, u.solLum)
            self.dl = WMAP9.luminosity_distance(self.z)

            #Rest frame GHz
            self.linename = ['12CO(3-2)','12CO(4-3)','12CO(5-4)','12CO(6-5)',
                             '12CO(7-6)','[CI](1-0)','[CI](2-1)']

            #W m^-2, observed
            self.linei0 = u.Quantity([1.2613e-20,2.19972e-20,4.76605e-20,
                                      5.01307e-20, 3.33625e-20,
                                      8.76716e-21,1.75036e-20], u.W / u.m**2)
            self.linei0 /= mu
            self.lineratio = numpy.ones_like(self.linei0)
            self.lineratio[-2:] = ciratio

        elif name == "GN20" :
            self.source = "Carilli et al. 2010"
            self.z = 4.055
            self.lir = u.Quantity(2.9e13, u.solLum)
            self.dl = WMAP9.luminosity_distance(self.z)

            self.linename = ['12CO(1-0)','12CO(2-1)','12CO(3-2)',
                             '12CO(4-3)','12CO(5-4)','12CO(6-5)',
                             '12CO(7-6)','12CO(8-7)','12CO(9-8)',
                             '[CI](1-0)']
            #[CI] (1-0) is assumed to be 1/3 of 12CO(4-3)
            self.linei0   = numpy.array([1.599e-22, 9.746e-22, 2.4e-21,
                                         4.57e-21,  8.374e-21, 8.22e-21,
                                         5.33e-21,  2.13e-21,  2.74e-21,
                                         1.52e-21])
            self.lineratio = numpy.ones_like(self.linei0)
            self.lineratio[-1:] = ciratio
        elif name == "FLS3" :
            mu = 1.25 # Mag limit
            self.source = "Riechers et al. 2013"
            self.z = 6.334
            self.lir = u.Quantity(4.16e13 / mu, u.solLum)
            self.dl = WMAP9.luminosity_distance(self.z)

            #First 6 lines are measured, rest are model (Bouthwell ratios,
            # or Arp220 beyond those), and a [CI](1-0) = 1/3 12CO(4-3) assump
            self.linename = ['12CO(1-0)','12CO(2-1)','12CO(6-5)',
                             '12CO(7-6)','12CO(9-8)','12CO(10-9)',
                             '12CO(3-2)','12CO(4-3)','12CO(5-4)',
                             '12CO(8-7)','[CI](1-0)']
            self.linei0   = u.Quantity([3.88e-23,  3.303e-22, 8.617e-21,
                                        1.064e-20, 1.305e-20, 2.049e-20,
                                        4.309e-21, 5.108e-21, 7.600e-21,
                                        1.143e-20, 1.703e-21], u.W / u.m**2)/mu
            self.lineratio = numpy.ones_like(self.linei0)
            self.lineratio[-1:] = ciratio

        else :
            raise KeyError("Unknown source %s" % name)

        #Rest frame GHz
        self.linefreq = numpy.array([linefreq[nm] for nm in self.linename])


    def obs_line_freq(self, redshift) :
        """Returns observed line frequency for all lines"""
        return self.linefreq / (1 + redshift)

    def linestrength(self, lir, redshift, lirslope=0.9, applyciratio=True) :
        """Predicts all available line strengths at (scalar) redshift 
        based on template and L_IR (in L sun if not specified)."""

        # Correction for luminosity distance
        dl_factor = (self.dl / WMAP9.luminosity_distance(redshift)).value**2

        # Determine normalization factor
        if isinstance(lir, u.Quantity):
            normfac = (lir.to(u.solLum) / self.lir).value**lirslope
        else:
            # Just assume lir is L_sun
            normfac = (lir / self.lir.value)**lirslope
            
        pred = dl_factor * normfac * self.linei0
        if applyciratio : 
            pred *= self.lineratio

        return (self.obs_line_freq(redshift), pred)

