import numpy as np
import astropy.units as u

"""Determines the sensitivity per bandwidth for ALMA
observations assuming the default conditions appropriate
for each band.

Per band, these are:
Band   Freq Range     Octile
----   ----------     -----------------
3      84-116         7th (5.186mm PWV)
4      125-163        6th (2.748mm PWV)
6      211-275        4th (1.262mm PWV)
7      275-373        4th (1.262mm PWV)
8      385-500        2nd (0.648mm PWV)
9      602-720        1st (0.472mm PWV)
"""

__all__ = ["band"]

band_id = np.array([3, 4, 6, 7, 8, 9], dtype=np.int16)
band_lowfreq = np.array([84.0, 125.0, 211.0, 275.0, 385.0, 602.0])
band_hifreq = np.array([116.0, 163.0, 275.0, 373.0, 500.0, 720.0])

class band:
    def __init__(self, freq):
        """ Object representing a single ALMA band.  Currently,
        only bands 3, 4, and 6 are supported.

        freq : float or astropy.units.Quantity
        ---------------------------------------
          Representative frequency.  Used to select the band.
          Assumed in GHz if not a Quantity.
        """

        if isinstance(freq, u.Quantity):
            if not freq.isscalar:
                raise ValueError("Only scalar frequency values supported")
            freq_ghz = freq.to(u.GHz).value
        else:
            freq_ghz = float(freq) # Assume GHz
        
        wband = np.nonzero((freq_ghz >= band_lowfreq) & 
                           (freq_ghz <= band_hifreq))[0]
        if len(wband) == 0:
            raise ValueError("Specified frequency is not in any known "\
                             "ALMA band")
        else:
            self._bandnum = band_id[wband[0]]

        self._minfreq = u.Quantity(band_lowfreq[wband[0]], u.GHz)
        self._maxfreq = u.Quantity(band_hifreq[wband[0]], u.GHz)
        self._loadtsys(self._bandnum)

    def _loadtsys(self, bandnum):
        """ Load the system temperature for a given band"""
        from pkg_resources import resource_filename
        import os.path
        from scipy.interpolate import interp1d

        # Construct tsys filename
        tsys_base = 'resources/band{:0d}_tsys.txt'.format(bandnum)
        tsys_filename = resource_filename(__name__, tsys_base)
        
        # Make sure it exists
        if not os.path.isfile(tsys_filename):
            raise IOError("Unable to find Tsys file "\
                          "for band {:d}".format(bandnum))
                        
        dat = np.loadtxt(tsys_filename)
        self._tsys_freq = dat[:, 0]
        self._tsys_tsys = dat[:, 1]
        self._min_tsys_freq = self._tsys_freq.min()
        self._max_tsys_freq = self._tsys_freq.max()
        self._tsys_interp = interp1d(self._tsys_freq, self._tsys_tsys,
                                     kind='linear')

    @property
    def bandnum(self):
        """ Gets ALMA band number"""
        return self._bandnum

    @property
    def minfreq(self):
        """ Minimum frequency of band"""
        return self._minfreq

    @property
    def maxfreq(self):
        """ Maximum frequency of band"""
        return self._maxfreq

    def tsys(self, freq):
        """ Get system temperature at specified frequency"""

        if isinstance(freq, u.Quantity):
            return u.Quantity(self._tsys_interp(freq.to(u.GHz).value), u.K)
        else:
            return u.Quantity(self._tsys_interp(freq), u.K)
        
    def aeff(self, freq) :
        """ Get effective area of 12m antennae at specified frequency"""

        # Note we don't check that we are within band
        prefac = -16 * np.pi**2
        #Surface rms (25 um) in GHz equivalent
        surffreq = 299792458.0 / 25e3
        #Effective area of 12ms
        if isinstance(freq, u.Quantity):
            f = freq.to(u.GHz).value
        else:
            f = freq

        return u.Quantity(113.1 * 0.72 *\
                          np.exp(prefac*(f / surffreq)**2), u.m**2)

    def sens(self, freq, tint = u.Quantity(1, u.s), 
             deltanu=u.Quantity(31.25, u.MHz), 
             nant=34, npol=2) :
        """ Computes 12m array sensitivity.
        
        Parameters
        ----------
        freq : astropy.units.Quantity
          Desired frequency for sensitivity.

        tint : astropy.units.Quantity
          Integration time.

        deltanu : astropy.units.Quantity
          Bandwidth over which sensitivity is defined.

        nant : int
          Number of 12m antennae

        npol : int
          Number of polarizations desired (so 2 means total power).

        Returns
        -------
        sens : astropy.units.Quantity
          Sensitivity. 

        Notes
        -----
          These compuations assume the OT default quantile for the
        band for a source at zenith.
        """

        import astropy.constants as const
        import math
    
        # Check bounds
        if isinstance(freq, u.Quantity):
            minf = freq.min()
            maxf = freq.max()
        elif isinstance(freq, np.ndarray):
            minf = u.Quantity(freq.min(), u.GHz)
            maxf = u.Quantity(freq.max(), u.GHz)
        else:
            minf = maxf = u.Quantity(float(freq), u.GHz)

        if minf < self._minfreq:
            raise ValueError("Frequency out of bounds (low) for band")
        if maxf > self._maxfreq:
            raise ValueError("Frequency out of bounds (high) for band")

        k_B = u.Quantity(const.k_B.si).value # J / K
        #Effective area of 12ms
        aeff = self.aeff(freq).value # m^2
        #Quantization efficiency
        eta_q = 0.96
        #Correlator efficiency
        eta_c = 0.88
        Tsys = self.tsys(freq).value
        # deltanu
        if isinstance(deltanu, u.Quantity):
            dnu_mhz = deltanu.to(u.MHz).value
        else:
            dnu_mhz = deltanu
        # tint
        if isinstance(tint, u.Quantity):
            tint_s = tint.to(u.s).value
        else:
            tint_s = tint
        
        return u.Quantity(1e26 * 2 * k_B * Tsys /\
                          (eta_q * eta_c * aeff * 
                           math.sqrt(nant*(nant-1) * npol * dnu_mhz * tint_s)),
                           u.mJy)

    def __repr__(self):
        return "ALMA band {0:d}".format(self.bandnum)
