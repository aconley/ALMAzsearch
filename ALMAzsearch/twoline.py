import numpy as np
import astropy.constants as cnst
import astropy.units as u

from .tuning import scan_tuning
from .line_pred import line_template

__all__ = ["twoline"]

""" Computes the S/N of the two most strongly detected lines given
an ALMA scan tuning"""

class twoline(object):
    def __init__(self, template_names, ciratio=1.0):
        if isinstance(template_names, str):
            # Just one
            self._ntemplates = 1
            self._templates = [line_template(template_names, ciratio=ciratio)]
            self._template_names = [template_names]
        else:
            self._ntemplates = len(template_names)
            self._templates = [line_template(nm, ciratio=ciratio) for
                               nm in template_names]
            self._template_names = template_names

    def sn(self, tuning, z, base_tint = u.Quantity(1, u.s),
           lir=u.Quantity(1e13, u.solLum),
           linewidth=u.Quantity(500, u.km / u.s), applyciratio=False):
        """Compute the S/N of the two most strongly detected lines
        for an input vector of z values as observed by tuning.
        If multiple templates are set, takes the one with the lowest
        S/N (over both lines).  Returns tuple of number of lines detected,
        S/N for each."""

        nz = len(z)
        # Set up return variables
        snret = np.zeros((nz, 2))
        nlines = np.zeros(nz)

        # Convert from velocity to frequency
        lw_factor = (linewidth / u.Quantity(cnst.c)).decompose().value

        for zidx in range(nz):
            curr_z = z[zidx]
            # Get predicted line strengths
            for i in range(self._ntemplates):
                tpl = self._templates[i]
                freq, i0 = tpl.linestrength(lir, curr_z, 
                                            applyciratio=applyciratio)
                # Figure out which lines are covered
                ncov = tuning.ncover(freq)
                # Only keep those that are covered
                nnonzero = np.count_nonzero(ncov)
                if nnonzero == 0:
                    # None observed
                    continue
                
                wnonzero = np.nonzero(ncov)[0]
                freq = freq[wnonzero]
                i0 = i0[wnonzero]
                
                # Get line fwhm in frequency
                linewidths = lw_factor * freq
                
                # Get predicted sensitivity for all lines
                # in W / m^2.  Note they come back in mJy,
                # so we have to multiply by the line width (in Hz)
                # to get W / m^2.  We estimate the sensitivity over
                # the linewidth provided by the user.
                # So, 1e-29 from mJy to W/m^2/Hz, then 1e9 because
                # the frequencies are in GHz, so a prefactor of 1e-20
                sens = np.array([1e-20 * lw.value *\
                                 tuning.sens(fq, base_tint=base_tint, 
                                             deltanu=lw).value
                                 for fq, lw in zip(freq, linewidths)])
                curr_sn = (i0 / sens).value
                
                if len(curr_sn) == 1:
                    # We only found one line.  This can only replace
                    # the previous best if that is also 0 or 1 line
                    sn0 = curr_sn[0]
                    if nlines[zidx] == 0: # New
                        nlines[zidx] = 1
                        snret[zidx, 0] = sn0
                    elif nlines[zidx] == 1 and sn0 > snret[zidx, 1]:
                        nlines[zidx] == 1
                        snret[zidx, 0] = sn0
                else:
                    # We found more than one line.  Take the two highest SN
                    # Base acceptance on lower of the two SN values
                    curr_sn.sort()
                    curr_sn1 = curr_sn[-1]
                    curr_sn2 = curr_sn[-2]
                    # If previous is 0 or 1 lines, always replace
                    # Otherwise accept if lower S/N is better than previous
                    if nlines[zidx] < 2 or curr_sn2 > snret[zidx, 2]:
                        snret[zidx, 0] = curr_sn1
                        snret[zidx, 1] = curr_sn2
                        nlines[zidx] = 2

        return (nlines, snret)
                        
    def __repr__(self):
        nm = ", ".join(self._template_names)
        return "ALMA two-line S/N calculator with templates: {:s}".format(nm)
