# Astropy unit equivalencies

import astropy.units as u

__all__ = ["radio_units"]

def radio_units(freq):
    """ Returns equivalencies for converting between si units (W m^-2)
    and radio units (Jy km s^-1).  Freq is the frequency of the observation."""

    if isinstance(freq, u.Quantity):
        if not freq.isscalar:
            raise ValueError("Freq must be scalar")
        freq_ghz = freq.to(u.GHz).value
    else:
        freq_ghz = float(freq)
        
    return [((u.Jy * u.km / u.s), (u.W / u.m**2), 
             lambda x: x * freq_ghz / 299792458e14,
             lambda x: 299792458e14 * x / freq_ghz)]
