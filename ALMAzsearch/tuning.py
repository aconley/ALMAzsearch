import numpy as np
import astropy.units as u

from .alma_sens import band

__all__ = ["scan_tuning"]

# This is the size of the 'gap' between the sideband sets
# don't support band 9 -- it's funny (DSB rather than 2SB)
# Note that band 6 has a funny split.
alma_split = {3: 8, 4: 8, 6: 12, 7: 8, 8: 8}


def overlap(x1, x2, y1, y2):
    """ Tests if [x1, x2] and [y1, y2] overlap, assuming
    x1 < x2, y1 < y2"""
    return (x1 <= y2) and (y1 <= x2)

###################################

class scan_tuning:
    """ Holds information about a single ALMA frequency scan tuning.

    Parameters
    ----------
    lofreq : astropy.units.Quantity
      The minimum frequency covered.  Note that this is not
      the SPW1 value in the OT (which is the center of the lowest tuning,
      so 1.875/2 GHz higher.

    base_tint : astropy.units.Quantity
      The integration time for a single tuning

    ntune : int
      The number of tunings.  Values of 1-5 are supported (since that's
      what the Cycle 2 ALMA OT allows).
    """

    def __init__(self, lofreq, base_tint, ntune=5, nantennae=34, npol=2):

        self._band = band(lofreq)
        self._nantennae = nantennae
        self._npol = npol

        if isinstance(base_tint, u.Quantity):
            if not base_tint.isscalar:
                raise ValueError("base_tint must be scalar")
            self._base_tint = base_tint.to(u.s)
        else:
            self._base_tint = u.Quantity(base_tint, u.s)

        if self._nantennae <= 0:
            raise ValueError("Invalid (non-positive) number of antennae")
        if not self._npol in {1, 2}:
            raise ValueError("Invalid number of polarizations")

        if isinstance(lofreq, u.Quantity):
            self.lowfreq = lofreq.to(u.GHz)
        else:
            self.lowfreq = u.Quantity(lofreq, u.GHz)

        if not self._band.bandnum in alma_split:
            errstr = "Band {0:d} not supported".format(self._band.bandnum)
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

        self._minfreq = self._min1.min()
        self._maxfreq = self._max2.max()

        # Make sure we are within range
        # Low edge already tested by band constructor
        if self._max2.max() > self._band.maxfreq:
            raise ValueError("Tuning extends beyond band boundary high")

    @property
    def ntune(self):
        return self._ntune
        
    @property
    def bandnum(self):
        return self._band.bandnum

    @property
    def nantennae(self):
        return self._nantennae

    @property
    def npol(self):
        return self._npol

    @property
    def base_tint(self):
        return self._base_tint

    @base_tint.setter
    def base_tint(self, value):
        if isinstance(value, u.Quantity):
            if not value.isscalar:
                raise ValueError("Integration time must be scalar")
            self._base_tint = value.to(u.s)
        else:
            self._base_tint = u.Quantity(value, u.s)
        

    def freq_coverage(self, idx):
        """ Returns the coverage for a scan index (range 0 to ntune-1)
        as min 1, max 1, min 2, max 2.  So, a frequency that is either
        between min 1 and max 1 or between min 2 and max 2 would be covered."""
        return u.Quantity([self._min1[idx].value, self._max1[idx].value,
                           self._min2[idx].value, self._max2[idx].value],
                          u.GHz)

    @property
    def min1(self):
        """ Minimum frequency of first sideband for each tuning"""
        return self._min1.copy()

    @property
    def max1(self):
        """ Maximum frequency of first sideband for each tuning"""
        return self._max1.copy()

    @property
    def min2(self):
        """ Minimum frequency of second sideband for each tuning"""
        return self._min2.copy()

    @property
    def max2(self):
        """ Maximum frequency of second sideband for each tuning"""
        return self._max2.copy()
        
    @property
    def minfreq(self):
        """ Minimum frequency covered; note coverage may not be contiguous"""
        return self._minfreq

    @property
    def maxfreq(self):
        """ Maximum frequency covered; note coverage may not be contiguous"""
        return self._maxfreq
        
    def overlap(self, tuning):
        """ Do any of the frequencies in this tuning overlap with those
        of another tuning?"""

        # Quick test for extremities
        if not overlap(self.minfreq, self.maxfreq, 
                       tuning.minfreq, tuning.maxfreq):
            return False

        # They may overlap -- we have to do the detailed comparison
        # of each individual tunings
        for idx1 in range(self._ntune):
            for idx2 in range(tuning._ntune):
                if overlap(self._min1[idx1], self._max1[idx1],
                           tuning._min1[idx2], tuning._max1[idx2]):
                    return True
                if overlap(self._min1[idx1], self._max1[idx1],
                           tuning._min2[idx2], tuning._max2[idx2]):
                    return True
                if overlap(self._min2[idx1], self._max2[idx1],
                           tuning._min1[idx2], tuning._max1[idx2]):
                    return True
                if overlap(self._min2[idx1], self._max2[idx1],
                           tuning._min2[idx2], tuning._max2[idx2]):
                    return True
        return False

    def ncover(self, freq):
        """ Find the number of tunings covering a specified frequency.

        Parameters
        ----------
        freq : astropy.units.Quantity
          Frequency.

        Returns
        -------
          The number of tunings covering the specified frequencies.
        """

        if isinstance(freq, u.Quantity):
            if not freq.isscalar:
                n1 = [np.count_nonzero((f >= self._min1) & (f <= self._max1))
                  for f in freq]
                n2 = [np.count_nonzero((f >= self._min2) & (f <= self._max2))
                  for f in freq]
                n1 = np.array(n1)
                n2 = np.array(n2)
            else:
                n1 = np.count_nonzero((freq >= self._min1) & 
                                      (freq <= self._max1))
                n2 = np.count_nonzero((freq >= self._min2) & 
                                      (freq <= self._max2))
        elif isinstance(freq, np.ndarray):
            min1, max1 = self._min1.value, self._max1.value
            min2, max2 = self._min2.value, self._max2.value
            if len(freq.shape) > 0:
                n1 = [np.count_nonzero((f >= min1) & (f <= max1)) for f in freq]
                n2 = [np.count_nonzero((f >= min2) & (f <= max2)) for f in freq]
                n1 = np.array(n1)
                n2 = np.array(n2)
            else:
                n1 = np.count_nonzero((freq >= min1) & (freq <= max1))
                n2 = np.count_nonzero((freq >= min2) & (freq <= max2))
        else:
            f = u.Quantity(float(freq), u.GHz)
            n1 = np.count_nonzero((f >= self._min1) & (f <= self._max1))
            n2 = np.count_nonzero((f >= self._min2) & (f <= self._max2))
        
        return n1 + n2

    def integration_time(self, freq):
        """ Return integration time at frequency"""
        return self._base_tint * self.ncover(freq)

    def sens(self, freq, deltanu=u.Quantity(31.25, u.MHz)):
        """ Get the sensitivity for a given frequency in the specified
        bandwidth (deltanu). Returns inf if the frequency is not covered by this
        setup."""
        
        if isinstance(freq, u.Quantity):
            if not freq.isscalar:
                raise ValueError("Only scalar frequency values supported")
            f = freq
        else:
            f = u.Quantity(float(freq), u.GHz)

        ncov = self.ncover(f)
        if ncov == 0:
            return u.Quantity(np.inf, u.mJy)

        tint = ncov * self._base_tint
        return self._band.sens(f, tint=tint, deltanu=deltanu, 
                               nant=self._nantennae,
                               npol=self._npol)
                
    def __repr__(self):
        rstr = "ALMA band {:d} frequency scan with {:d} tunings, "\
               "min/max frequency: {:s}/{:s} and base integration: {:s}"
        return rstr.format(self.bandnum, self._ntune, self._minfreq, 
                           self.maxfreq, self.base_tint)
