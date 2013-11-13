import numpy as np
import astropy.constants as cnst
import astropy.units as u

from .tuning import scan_tuning
from .line_pred import line_template

__all__ = ["twoline"]

""" Computes the S/N of the two most strongly detected lines given
an ALMA scan tuning"""

class twoline(object):
    """ An object which can compute S/N ratios provided with an ALMA 
    observation

    Parameters
    ----------
    template_names : str or [str]
      Names of line-strength templates to use in the calculation

    ciratio : float
      Bonus multiplier to returned S/N of [CI] lines.  This is to
      account for the fact that they occur in a fixed relation
      with nearby 12CO lines, and hence a lower S/N is acceptable when
      searching for them.

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

    def __init__(self, template_names, ciratio=1.0, Td=u.Quantity(55, u.K), 
                 beta=1.8, lambda0=u.Quantity(200.0, u.um), alpha=4):
        if isinstance(template_names, str):
            # Just one
            self._ntemplates = 1
            self._templates = [line_template(template_names, ciratio=ciratio,
                                             Td=Td, beta=beta, lambda0=lambda0,
                                             alpha=alpha)]
            self._template_names = [template_names]
        else:
            self._ntemplates = len(template_names)
            self._templates = [line_template(nm, ciratio=ciratio, Td=Td, 
                                             beta=beta, lambda0=lambda0,
                                             alpha=alpha) for
                               nm in template_names]
            self._template_names = template_names

    def _get_sn(self, tuning, freq, i0, lw_factor):
        # Figure out which lines are covered
        ncov = tuning.ncover(freq)

        # Only keep those that are covered
        nnonzero = np.count_nonzero(ncov)
        if nnonzero == 0:
            return None
                
        wnonzero = np.nonzero(ncov)[0]
        wfreq = freq[wnonzero]
        wi0 = i0[wnonzero]
                
        # Get line fwhm in frequency
        linewidths = lw_factor * wfreq
                
        # Get predicted sensitivity for all lines
        # in W / m^2.  Note they come back in mJy,
        # so we have to multiply by the line width (in Hz)
        # to get W / m^2.  We estimate the sensitivity over
        # the linewidth provided by the user.
        # So, 1e-29 from mJy to W/m^2/Hz, then 1e9 because
        # the frequencies are in GHz, so a prefactor of 1e-20
        sens = np.array([1e-20 * lw.value * tuning.sens(fq, deltanu=lw).value
                         for fq, lw in zip(wfreq, linewidths)])
        return (wi0 / sens).value
                
    def sn(self, tunings, z, norm=u.Quantity(3e13, u.solLum),
           linewidth=u.Quantity(500, u.km / u.s), lirslope=0.9,
           obs_wave=u.Quantity(500, u.um), applyciratio=False, minsn=None):
        """ Compute the S/N values.

        Parameters
        ----------
        tunings : scan_tuning or [scan_tuning]
          An object or set of objects specifying the ALMA tunings to use
          in the computation.  These must not overlap -- any overlaps are
          ignored.

        z : float or np.ndarray
          Redshifts to compute S/N values at

        norm : astropy.units.Quantity
          Quantity to use for normalization purposes.  If in power
          units (e.g., L_sun) then assumed to be L_IR (8-1000um).  
          If in spectral flux density units (e.g., mJy) then 
          assumed to be observed continuum flux density at obs_wave.

        linewidth : astropy.units.Quantity
          Assumed linewidth in km/s.

        lirslope : float
          Log slope of L_IR, L_line relation, such that 
          L_line \propto L_IR^lirslope

       obs_wave : astropy.units.Quantity
          If using flux density based normalization, observer frame
          wavelength of observation.

        applyciratio : bool
           Apply [CI] S/N bonus in computation.

        minsn : float or None
           Require minimum S/N before recording.

        Returns
        -------
        A tuple of (nlines, sn) where nlines is a integer np.ndarray of length
        nz (=number of z values) giving the number of lines covered by
        the frequency setup, and sn is a (nz, 2) np.ndarray giving
        the recovered S/N of the two highest S/N lines, such that [:, 1]
        is the lower of the two S/N values.  If only one line is detected,
        [:,1] is 0, but [:,0] is filled.
        
        Notes
        -----
        At each redshift, for each template, computes the S/N of the two 
        most strongly detected lines covered by tunings.  Then takes the
        lowest S/N across that set of templates (e.g., the pessimistic case)
        for each redshift.

        See the documentation for `line_template.linestrength` for discussion 
        of the normalization.
        """

        # Test if any tunings overlap
        try:
            ntune = len(tunings)
            if ntune > 1:
                for idx1 in range(ntune-1):
                    for idx2 in range(idx1 + 1, ntune):
                        if tunings[idx1].overlap(tunings[idx2]):
                            errstr = "Tunings {:d} and {:d} overlap"
                            raise Exception(errstr.format(idx1, idx2))
        except TypeError:
            # This just means there was only one
            pass

        nz = len(z)
        # Set up return variables
        snret = np.zeros((nz, 2))
        nlines = np.zeros(nz, dtype=np.int16)

        # Convert from velocity to frequency
        lw_factor = (linewidth / u.Quantity(cnst.c)).decompose().value

        for zidx in range(nz):
            curr_z = z[zidx]
            # Get predicted line strengths
            for i in range(self._ntemplates):
                tpl = self._templates[i]
                lnname, freq, i0 = tpl.linestrength(norm, curr_z,
                                                    lirslope=lirslope,
                                                    obs_wave=obs_wave)
                if applyciratio:
                    i0 *= tpl.lineratio

                # For each tuning, get the S/N values
                try:
                    curr_sn = [self._get_sn(tune, freq, i0, lw_factor)
                               for tune in tunings]
                    # We have to flatten.  First pick off the Nones
                    curr_sn = list(filter(lambda x: not x is None,
                                          curr_sn))
                    if len(curr_sn) == 0:
                        continue
                    curr_sn = np.concatenate(curr_sn)
                except TypeError:
                    curr_sn = self._get_sn(tunings, freq, i0, lw_factor)
                    if curr_sn is None:
                        continue

                # Process these out.
                nl = len(curr_sn)
                if nl == 1:
                    # We only found one line.  If we already have a two-line
                    # solution at this redshift (from a different template),
                    # ignore this.  But if we have a one-line, take this
                    # if it is -lower- S/N (remember, we are being pessimistic)
                    # Of course, if we have a no-line solution, always take
                    sn0 = curr_sn[0]
                    if not minsn is None:
                        if sn0 < minsn:
                            continue
                    if nlines[zidx] == 0: # New
                        nlines[zidx] = 1
                        snret[zidx, 0] = sn0
                    elif nlines[zidx] == 1 and sn0 < snret[zidx, 0]:
                        snret[zidx, 0] = sn0
                else:
                    # We found more than one line.  Take the two highest SN
                    # Base acceptance on lower of the two SN values
                    curr_sn.sort()
                    curr_sn1 = curr_sn[-1]
                    curr_sn2 = curr_sn[-2]
                    if not minsn is None:
                        if curr_sn2 < minsn:
                            continue
                    # If previous is 0 or 1 lines, always replace
                    # Otherwise accept if lower S/N is -worse- than previous
                    if nlines[zidx] < 2 or curr_sn2 < snret[zidx, 1]:
                        snret[zidx, 0] = curr_sn1
                        snret[zidx, 1] = curr_sn2
                        nlines[zidx] = nl

        return (nlines, snret)
                        
    def __repr__(self):
        nm = ", ".join(self._template_names)
        return "ALMA two-line S/N calculator with templates: {:s}".format(nm)
