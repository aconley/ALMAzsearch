import numpy as np
import astropy.units as u
from astropy.cosmology import WMAP9

"""Predicts line intensities at a given redshift from a template"""

from .co_lines import linefreq
from .modified_blackbody import modified_blackbody

__all__ = ["line_template"]

class line_template(object) :
    """ Object for predicting observed line strengths
    
    Two normalization modes are available: L_IR (8-1000) based,
    and observed flux density based.  The latter is based on an
    optically thick modified blackbody model.

    Parameters
    ----------
    name : str
      Name of the template.  The currently supported values are
      Arp220, Eyelash, HLSW01, FLS3, GN20, and ID141

    ciratio : float
      S/N bonus factor for [CI] lines -- see twoline documentation.

    Td : astropy.units.Quantity
      Dust temperature (rest frame) for observed flux normalization.

    beta : float
      Dust beta for observed flux normalization.

    lambda0 : astropy.units.Quantity
      Rest frame wavelength where dust becomes optically thick for
      flux normalization.

    alpha : float
      Wein-side power-law slope for flux normalization.
    """


    def __init__(self, name, ciratio=1.0, Td=u.Quantity(55, u.K), 
                 beta=1.8, lambda0=u.Quantity(200.0, u.um), alpha=4):
        # We store which lines are observed (or predicted) as
        # well as L_IR and the observed line strengths (both
        # corrected for lensing)

        # Most templates are combinations of observed lines plus model
        # predictions.  For objects without [CI](1-0), a value of
        # I_[CI](1-0) = 0.39 I_12CO(4-3), which is a reasonably representative
        # value for the things in our sample with actual measurements.
        # The two exceptions are the Cosmic Eyelash (with very strong [CI]
        # lines) and the SPT composite, which has I_[CI](1-0) = 1/3 I_12CO(4-3)
        # For objects without [CI](2-1), it is assumed that 
        # I_[CI](2-1) = 5 / 3 I_[CI](1-0)

        if name == "Arp220" :
            self._source = "Rangwala et al. (2011)"
            self._z = 0.0181
            self._dl = u.Quantity(77, u.Mpc)
            self._lir = u.Quantity(1.77e12, u.solLum)
            
            #Names of lines
            self._linename = ['12CO(1-0)','12CO(2-1)','12CO(3-2)',
                              '12CO(4-3)','12CO(5-4)','12CO(6-5)',
                              '12CO(7-6)','12CO(8-7)','12CO(9-8)',
                              '12CO(11-10)','12CO(12-11)','12CO(13-12)',
                              '[CI](1-0)','[CI](2-1)']
            
            
            #W m^-2, observed
            self._linei0 = u.Quantity([1.582e-18,  8.512e-18,  4.162e-17,
                                       6.9895e-17, 7.0314e-17, 9.3966e-17,
                                       9.3205e-17, 10.014e-17, 10.1133e-17,
                                       5.313e-17,  3.803e-17,  2.207e-17,
                                       3.015e-17, 5.73e-17], u.W / u.m**2)
            

            self.lineratio = np.ones_like(self._linei0.value)
            self.lineratio[-2:] = ciratio

        elif name == "Eyelash" or name == "SMMJ2135-0102":
            self._source = "Danielson et al. (2011)"
            self._z = 2.3259
            self._dl = WMAP9.luminosity_distance(self._z)
            self._lir = u.Quantity(2.3e12, u.solLum)

            #Names of lines, all observed
            self._linename = ['12CO(1-0)','12CO(3-2)','12CO(4-3)',
                              '12CO(5-4)','12CO(6-5)','12CO(7-6)',
                              '12CO(8-7)','12CO(9-8)',
                              '[CI](1-0)','[CI](2-1)']
            
            #W m^-2, observed, corrected for lensing
            self._linei0 = u.Quantity([7.68e-23, 1.409e-21, 2.461e-21, 
                                       3.323e-21, 4.58e-21, 3.13e-21, 
                                       2.5e-21, 1.15e-21, 
                                       2.43e-21, 4.044e-21], 
                                      u.W / u.m**2)
            self.lineratio = np.ones_like(self._linei0.value)
            self.lineratio[-2:] = ciratio

        elif name == "HLSW01" or name == "LSW01" or name == "LSW_01" :
            self._source = "K.S. Scott et al. (2011)"
            self._z = 2.9575
            self._lir = u.Quantity(1.43e13, u.solLum) # Magnification corrected
            self._dl = WMAP9.luminosity_distance(self._z)

            # 2-1, 4-3, 6-5, 8-7 are RADEX model based, [CI](1-0) and 
            # [CI](2-1) as described above. Values are magnification corrected
            self._linename = ['12CO(1-0)', '12CO(2-1)', '12CO(3-2)',
                              '12CO(4-3)', '12CO(5-4)', '12CO(6-5)',
                              '12CO(7-6)', '12CO(8-7)', '12CO(9-8)',
                              '12CO(10-9)' ,'[CI](1-0)', '[CI](2-1)']
            self._linei0 = u.Quantity([1e-22,   1.0e-21, 2.6e-21,
                                       0.6e-20, 1.1e-20, 1.5e-20,
                                       2.2e-20, 1.7e-20, 1e-20, 
                                       1.3e-20, 2.3e-21, 3.9e-21],
                                      u.W / u.m**2)
            self.lineratio = np.ones_like(self._linei0.value)
            self.lineratio[-2:] = ciratio

        elif name == "ID141" :
            #Observer frame SED params; note none of the below are
            # corrected for lensing.  The CO (6-5) and (3-2) are from
            # their model.  Since the lensing magnification is not known,
            # a value of mu = 5 has been assumed
            mu = 5.0
            self._source = "Cox et al. (2011)"
            self._z = 4.243
            self._lir = u.Quantity(8.5e13 / mu, u.solLum)
            self._dl = WMAP9.luminosity_distance(self._z)

            #Rest frame GHz
            self._linename = ['12CO(3-2)', '12CO(4-3)', '12CO(5-4)',
                              '12CO(6-5)', '12CO(7-6)',
                              '[CI](1-0)', '[CI](2-1)']

            #W m^-2, observed
            self._linei0 = u.Quantity([1.2613e-20, 2.19972e-20, 4.76605e-20,
                                       5.01307e-20, 3.33625e-20,
                                       8.76716e-21,1.75036e-20], 
                                      u.W / u.m**2) / mu
            self.lineratio = np.ones_like(self._linei0.value)
            self.lineratio[-2:] = ciratio

        elif name == "GN20" :
            self._source = "Carilli et al. (2010)"
            self._z = 4.055
            self._lir = u.Quantity(2.9e13, u.solLum)
            self._dl = WMAP9.luminosity_distance(self._z)

            # 1-0, 2-1, 4-3, 5-4, 6-5 are observed
            # Other CO are from a two component LVG model, plus the
            # [CI] assumptions discussed above
            self._linename = ['12CO(1-0)', '12CO(2-1)', '12CO(3-2)',
                              '12CO(4-3)', '12CO(5-4)', '12CO(6-5)',
                              '12CO(7-6)', '12CO(8-7)', '12CO(9-8)',
                              '[CI](1-0)', '[CI](2-1)']
            self._linei0   = u.Quantity([1.60e-22, 9.75e-22, 2.74e-21,
                                         4.57e-21, 8.37e-21, 8.22e-21,
                                         5.33e-21, 2.13e-21, 2.74e-21,
                                         1.78e-21, 2.97e-21], u.W / u.m**2)
            self.lineratio = np.ones_like(self._linei0.value)
            self.lineratio[-2:] = ciratio
        elif name == "FLS3" or name == "HFLS3" or name == "FLS_3" :
            mu = 1.25 # Mag limit
            self._source = "Riechers et al. (2013)"
            self._z = 6.334
            self._lir = u.Quantity(4.16e13 / mu, u.solLum)
            self._dl = WMAP9.luminosity_distance(self._z)

            #First 7 lines are measured, rest of the CO lines (4-3, 5-4, 8-7)
            # are model (the RADEX model from the paper), and 
            # [CI] assumptions discussed above.
            self._linename = ['12CO(1-0)', '12CO(2-1)', '12CO(3-2)',
                              '12CO(6-5)', '12CO(7-6)','12CO(9-8)',
                              '12CO(10-9)', '12CO(4-3)','12CO(5-4)',
                              '12CO(8-7)','[CI](1-0)','[CI](2-1)']
            self._linei0   = u.Quantity([3.88e-23, 3.30e-22, 1.13e-21,
                                         8.62e-21, 1.06e-20, 1.31e-20, 
                                         2.05e-20, 2.73e-21, 4.85e-21,
                                         1.26e-20, 1.06e-21, 1.77e-21],
                                        u.W / u.m**2) / mu
            self.lineratio = np.ones_like(self._linei0.value)
            self.lineratio[-2:] = ciratio
        elif name == "SPTComposite":
            self._source = "Spilker et al. (2013)"
            self._z = 3.0
            self._lir = u.Quantity(5e13, u.solLum)
            self._dl = WMAP9.luminosity_distance(self._z)

            # 12CO(1-0) to 12CO(6-5) are observed. 12CO(7-6) and 12CO(8-7)
            # are based on the Spilker hot, diffuse model (fig 5).
            # [CI](1-0) is observed, [CI](2-1) using assumption discussed
            # above
            self._linename = ['12CO(1-0)', '12CO(2-1)', '12CO(3-2)',
                              '12CO(4-3)', '12CO(5-4)', '12CO(6-5)',
                              '12CO(7-6)', '12CO(8-7)', '[CI](1-0)',
                              '[CI](2-1)']
            self._linei0 = u.Quantity([6.92e-22, 6.10e-21, 1.61e-20,
                                       2.66e-20, 5.79e-20, 6.83e-20,
                                       5.10e-20, 4.26e-20, 8.57e-21,
                                       1.43e-20], u.W / u.m**2)
            self.lineratio = np.ones_like(self._linei0.value)
            self.lineratio[-2:] = ciratio
        elif name == "J16359+6612":
            # This is just the B image
            # This object has -extremely- poorly constrainted L_IR, so
            # is probably not suitable to actually use.
            mu = 22.0
            self._source = "Kneib et al. (2005), Weiss et al. (2005), Walter et al. (2011)"
            self._z = 2.5174
            # L_FIR to L_IR, mag corr already.  Very uncertain.
            self._lir = u.Quantity(1.46e12, u.solLum) * 1.3 
            self._dl = WMAP9.luminosity_distance(self._z)

            # All observed
            self._linename = ['12CO(3-2)','12CO(4-3)','12CO(5-4)',
                              '12CO(6-5)','12CO(7-6)','[CI](1-0)',
                              '[CI](2-1)']
            self._linei0 = u.Quantity([9.18e-21, 1.75e-20, 2.79e-20,
                                       2.62e-20, 1.91e-20, 7.93e-21,
                                       1.69e-20], u.W / u.m**2) / mu
            self.lineratio = np.ones_like(self._linei0.value)
            self.lineratio[-2:] = ciratio
            
        else:
            raise KeyError("Unknown source %s" % name)

        #Rest frame GHz
        self._linefreq = u.Quantity([linefreq[nm] for nm in self._linename],
                                   u.GHz)
        
                                     
        self._name = name

        # Set up dust model
        self._Td = Td.to(u.K).value # In kelvin
        self._beta = float(beta)
        self._lambda0 = lambda0.to(u.um).value # In microns
        self._alpha = float(alpha)
        # Get reference flux density for object
        # The idea is to compute the L_IR in the case this had
        # a 500um flux density of 1.0 mJy, then scale up to match
        # the desired L_IR and store the 500um flux density, recording
        # that 500um flux for reference.  So we are assuming that
        # continuum luminosities scale as L_IR
        opz = 1.0 + self._z
        sedmodel = modified_blackbody(self._Td / opz, self._beta,
                                      self._lambda0 * opz, self._alpha,
                                      1.0)
        lirprefac = 3.11749657e16 * self.dl.value**2  # To L_sun
        pred_lir = lirprefac * sedmodel.freq_integrate(8.0 * opz, 1000.0 * opz)
        self._f500 = self.lir.value / pred_lir # In mJy
        self._basesed = modified_blackbody(self._Td / opz, self._beta,
                                           self._lambda0 * opz, self._alpha,
                                           self._f500)

    @property
    def name(self):
        return self._name

    @property
    def source(self):
        return self._source

    @property
    def redshift(self):
        return self._z

    @property
    def lir(self):
        return self._lir

    @property
    def dl(self):
        return self._dl

    @property 
    def lines(self):
        return (self._linename, self._linefreq, self._linei0)

    def base_continuum(self, redshift, wave=u.Quantity(500, u.um)):
        """ Predicts the continuum flux for the template at a given redshift
        and observer frame wavelength.

        Parameters
        ----------
        redshift : float
          Redshift to compute expected brightness at.
        
        wave : u.Quantity
          Observer frame wavelength.
        """

        this_dl = WMAP9.luminosity_distance(redshift).value
        this_opz = 1.0 + redshift
        tpl_opz = 1.0 + self._z
        this_wave = wave.to(u.um).value 
        snu = this_opz / tpl_opz * (self._dl.value/this_dl)**2 * \
              self._basesed(this_wave * tpl_opz / this_opz)
        return u.Quantity(snu, u.mJy)

    def _lirfac(self, redshift, snu, wave):
        return snu.to(u.mJy).value / self.base_continuum(redshift, wave).value

    def lir_snu(self, redshift, snu, wave=u.Quantity(500, u.um)):
        """ Predict L_IR required in order to match observed snu at wave and
        redshift"""

        return self._lir * self._lirfac(redshift, snu, wave)

    def obs_line_freq(self, redshift) :
        """Returns observed line frequency for all lines"""
        return self._linefreq / (1 + redshift)

    def linestrength(self, norm, redshift, lirslope=0.9,
                     obs_wave=u.Quantity(500, u.um)):
        """ Predict line strengths
        
        Parameters
        ----------
        norm : astropy.units.Quantity
          Quantity to use for normalization purposes.  If in power
          units (e.g., L_sun) then assumed to be L_IR (8-1000um).  
          If in spectral flux density units (e.g., mJy) then 
          assumed to be observed continuum flux density at obs_wave.

        redshift : float
          Redshift to evaluate line strengths at.

        lirslope : float
          Log slope of L_IR, L_line relation, such that 
          L_line \propto L_IR^lirslope
        
        obs_wave : astropy.units.Quantity
          If using flux density based normalization, observer frame
          wavelength of observation.

        Returns
        -------
          A tuple of line name, observer frame frequency, and
          line strength.  The first is a list of strings, the others
          are astropy.units.Quantity objects.

        Notes
        -----
          The ciratio S/N bonus is not applied here.
          If the normalization is done in power units, then the rest frame line
          strengths are simply scaled up using lirslope.  If observed continuum
          flux density is used, the L_IR correction is computed so that
          the template continuum matches the observed value, then the
          same line strength adjustment is applied.
        """

        # Correction for luminosity distance
        dl_factor = (self._dl / WMAP9.luminosity_distance(redshift)).value**2

        # Determine normalization factor
        if isinstance(norm, u.Quantity):
            ptype = norm.unit.physical_type
            if ptype == "power":
                lirfac = norm.to(u.solLum).value / self._lir.value
            elif ptype == "spectral flux density":
                lirfac = self._lirfac(redshift, norm, obs_wave)
            else:
                raise ValueException("Normalization not power or flux density")
        else:
            # Just assume lir is in L_sun
            lirfac = float(norm) / self._lir.value

        normfac = lirfac ** lirslope
        pred = dl_factor * normfac * self._linei0
        return (self._linename, self.obs_line_freq(redshift), pred)

    def __repr__(self):
        rpr = "Line strength template for {:s} from {:s}"
        return rpr.format(self._name, self._source)
