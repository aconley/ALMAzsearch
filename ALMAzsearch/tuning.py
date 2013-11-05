import numpy as np
import astropy.units as u

from .alma_sens import band

__all__ = ["scan_tuning"]

# This is the size of the 'gap' between the sideband sets
# don't support band 9 -- it's funny (DSB rather than 2SB)
# Note that band 6 has a funny split.
alma_split = {3: 8, 4: 8, 6: 12, 7: 8, 8: 8}

class scan_tuning:
    """ Holds information about a single ALMA frequency scan tuning.

    lofreq is the minimum frequency covered.  Note that this is not
      the SPW1 value in the OT (which is the center of the lowest tuning,
      so 1.875/2 GHz higher.

    ntune is the number of tunings.  Values of 1-5 are supported (although
    1 would be quite silly).
    """

    def __init__(self, lofreq, ntune=5):
        self.band = band(lofreq)

        if isinstance(lofreq, u.Quantity):
            self.lowfreq = lofreq.to(u.GHz)
        else:
            self.lowfreq = u.Quantity(lofreq, u.GHz)

        if not self.band.bandnum in alma_split:
            errstr = "Band {0:d} not supported".format(self.band.bandnum)
            raise ValueError(errstr)

        # Make an array for each tuning giving the ranges covered
        if ntune <= 0:
            raise ValueError("Invalid (non-positive) number of tunings")
        if ntune > 5:
            raise ValueError("ALMA supports 5 or fewer tunings in spectral scan")
        self._ntune = ntune

        # This is the way that the OT sets up spectral scan mode.
        bwidth = 1.875 # GHz, with edges not used
        self._min1 = u.Quantity(self.lowfreq.value + 
                                2 * bwidth * np.arange(ntune), u.GHz)
        self._max1 = self._min1 + 2 * u.Quantity(bwidth, u.GHz)
        # This is the gap plus the sideband above the minimum frequency
        self._min2 = self._min1 + \
                     u.Quantity(alma_split[self.bandnum] + 4, u.GHz)
        self._max2 = self._min2 + 2 * u.Quantity(bwidth, u.GHz)

        # Make sure we are within range
        # Low edge already tested by band constructor
        if self._max2.max() > self.band.maxfreq:
            raise ValueError("Tuning extends beyond band boundary high")

    @property
    def ntune(self):
        return self._ntune
        
    @property
    def bandnum(self):
        return self.band.bandnum

    def freq_coverage(self, idx):
        """ Returns the coverage for a scan index (range 0 to ntune-1)
        as min 1, max 1, min 2, max 2.  So, a frequency that is either
        between min 1 and max 1 or between min 2 and max 2 would be covered."""
        return u.Quantity([self._min1[idx].value, self._max1[idx].value,
                           self._min2[idx].value, self._max2[idx].value],
                          u.GHz)

    @property
    def min1(self):
        return self._min1.copy()

    @property
    def max1(self):
        return self._max1.copy()

    @property
    def min2(self):
        return self._min2.copy()

    @property
    def max2(self):
        return self._max2.copy()

    def ncover(self, freq):
        """ Return the number of tunings covering a specified frequency"""

        if isinstance(freq, u.Quantity):
            if not freq.isscalar:
                raise ValueError("Only scalar frequencies supported")
            n1 = np.count_nonzero((freq >= self._min1) & (freq <= self._max1))
            n2 = np.count_nonzero((freq >= self._min2) & (freq <= self._max2))
        else:
            f = u.Quantity(float(freq), u.GHz)
            n1 = np.count_nonzero((f >= self._min1) & (f <= self._max1))
            n2 = np.count_nonzero((f >= self._min2) & (f <= self._max2))
        
        return n1 + n2

    def sens(self, freq, base_tint = u.Quantity(1, u.s),
             deltanu = u.Quantity(31.25, u.MHz),
             nant=34, npol=2) :
        """ Get the sensitivity for a given frequency in the specified
        bandwidth (deltanu) with the base integration time (for one tuning)
        base_tint.  Returns inf if the frequency is not covered by this
        setup."""
        
        if isinstance(freq, u.Quantity):
            if not freq.isscalar:
                raise ValueError("Only scalar frequency values supported")
            f = freq
        else:
            f = u.Quantity(float(freq), u.GHz)

        if isinstance(base_tint, u.Quantity):
            if not base_tint.isscalar:
                raise ValueError("Only scalar integration times supported")
            t = base_tint
        else:
            t = u.Quantity(float(base_tint), u.s)
            
        ncov = self.ncover(f)
        if ncov == 0:
            return u.Quantity(np.inf, u.mJy)

        tint = ncov * t
        return self.band.sens(f, tint=t, deltanu=deltanu, nant=nant,
                              npol=npol)
                
