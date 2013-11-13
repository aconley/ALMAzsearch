# Astropy unit equivalencies for radio line strength units

from math import pi
import astropy.units as u
import astropy.constants.si as _si

__all__ = ["radio_lines_simple", "radio_lines"]

def radio_lines_simple(disp_obs):
    """ Returns equivalence pairs between Jy km/s, a line-strength
    unit often reported in radio astronomy, and SI units.

    Parameters
    ----------
    disp_obs : `Quantity` with spectral units
        The `spectral` equivalent `Unit` (e.g., frequency or
        wavelength) in the observer frame.

    Notes
    -----
      This is a less-capable version of `radio_lines`, but does not
      require a redshift or cosmology.
    """
    
    if not disp_obs.isscalar:
        raise ValueError("disp_obs must be scalar")
    fobs = disp_obs.to(u.Hz, equivalencies=spectral()).value
    c_si = _si.c.value # m / s
    return [(u.Jy * u.km / u.s, u.W / u.m**2,
             lambda x: x * 1e-23 * fobs / c_si,
             lambda x: c_si * x * 1e23 / fobs)]
    

def radio_lines(disp_obs, disp_rest, lumdist):
    """
    Returns a list of equivalence pairs between observational
    and non-observational line strength units in astronomy.

    Observationally, line units are typically reported in Jy km/s.
    These are often converted to either solar luminosities or
    K km/s pc^2 by assuming a cosmological model.  These equivalencies
    allow the user to go between these pairs.

    Parameters
    ----------
    disp_obs : `Quantity` with spectral units
        The `spectral` equivalent `Unit` (e.g., frequency or
        wavelength) in the observer frame.

    disp_rest : `Quantity` with spectral units
        The `spectral` equivalent `Unit` (e.g., frequency or
        wavelength) in the rest frame.

    lumdist : `Quantity` with units of distance
        The luminosity distance to the source

    Notes
    -----
      This is a more capable version of `radio_lines_simple`, but at
      the price of requiring additional information.  If you
      just want to convert between W/m^2 and Jy km/s, you can use
      `radio_lines_simple`.
    """

    # The conversion from Jy km/s to W/m^2 is
    #     1 [Jy km/s] = 1e-26 * nu_obs / c_kms [W/m2]
    #  where 1e-26 is because of the definition of a jansky and c_kms
    #  is the speed of light in km/s

    # The conversion from W/m^2 to L_sun is
    #     1 [W/m^2] = 4 * pi * D_L^2 / Lsol [L_sun]
    # where D_L is the luminosity distance in m, and Lsol is the luminosity
    # of the sun in SI units

    # The conversion from K km/s pc^2 to W/m^2 is
    #     2 * 1000 * nu_obs**3 * k_B / (1+z) c^3 (D_A)**2
    # Where 2 nu^2/c^2 k_B comes from the R-J expression for a Blackbody
    # 1000 nu / c comes from km/s to frequency width
    # D_A is measured in pc because the expression is (L/D_A)**2 where
    #  L is the size of the emitting region, which is 1pc for K km/s pc^2

    # All the inter-conversoins follow from those.

    # Approximate versions (less precision on constants in front) 
    #  are given in Carilli and Walter, AARA, 2013, 51, section 2.4
    # Unfortunately, they don't derive the 'magic numbers' in front,
    # and I (@aconley) am not aware of anywhere where they are explained,
    # but they can be worked out by hand in terms of fundamental constants 
    # and unit conversions.

    import astropy.cosmology

    if not disp_obs.isscalar:
        raise ValueError("disp_obs must be scalar")
    fobs = disp_obs.to(u.Hz, equivalencies=spectral()).value
    if not disp_rest.isscalar:
        raise ValueError("disp_rest must be scalar")
    frest = disp_obs.to(u.Hz, equivalencies=spectral()).value
    opz = frest / fobs

    dl = lumdist.to(astrophys.Mpc)
    if not dl.isscalar:
        raise ValueError("User provided luminosity distance must be scalar")
    if dl.value <= 0.0:
        raise ValueError("Invalid (non-positive) luminosity distance")

    # Set up quantities we will reuse below
    dl_mpc = dl.value # Mpc, convenient to keep around because of the pc^2
    dl_m = dl.to(u.m).value # m
    # SI 4pi D_L^2 / L_sun 
    fpidl2overLsun_si = 4.0 * pi * dl_m ** 2 / _si.L_sun.value 
    c_si = _si.c.value # m/s
    c2overkb = c_si ** 2 / _si.k_B.value # c^2 / Boltzmann in SI units
    c3overkb = c_si * c2overkb # c^3 / Boltzmann in SI units

    # This is L_sun * c**3 / (4 pi (1pc/1m)**2 * k_B) in SI units,
    # which shows up in the lsun to K km/s pc^2 conversion
    lsun_kkmspc2_const = _si.L_sun.value * c3overkb /\
                         (4.0 * pi * astrophys.parsec.to(u.m) ** 2)

    def jykms_to_wm2(x):
        return x * 1e-23 * fobs / c_si
    
    def wm2_to_jykms(x):
        return c_si * x * 1e23 / fobs

    def wm2_to_lsun(x):
        return x * fpidl2overLsun_si

    def lsun_to_wm2(x):
        return x / fpidl2overLsun_si

    def jykms_to_lsun(x):
        return x * 1e-23 * fobs * fpidl2overLsun_si / c_si

    def lsun_to_jykms(x):
        return 1e23 * x * c_si / (fobs * fpidl2overLsun_si)

    def wm2_to_kkmspc2(x):
        return 0.5e9 * x * c3overkb * dl_mpc**2 / (frest ** 3)
        
    def kkmspc2_to_wm2(x):
        return 2e-9 * x * frest ** 3 / (c3overkb * dl_mpc ** 2)

    def jykms_to_kkmspc2(x):
        return 0.5e-14 * x * c2overkb * dl_mpc**2 / (opz * frest ** 2)

    def kkmspc2_to_jykms(x):
        return 2e14 * x * opz * frest**2 / (c2overkb * dl_mpc**2)

    def lsun_to_kkmspc2(x):
        return x * 5e-4 * lsun_kkmspc2_const * frest ** (-3)

    def kkmspc2_to_lsun(x):
        return x * 2e3 * frest ** 3 / lsun_kkmspc2_const 
    
    wm2 = u.W / u.m**2
    jykms = astrophys.Jy * u.km / u.s
    kkmspc2 = u.K * u.km * u.pc**2 / u.s
    
    return [
        (wm2, jykms, wm2_to_jykms, jykms_to_wm2),
        (wm2, astrophys.solLum, wm2_to_lsun, lsun_to_wm2),
        (wm2, kkmspc2, wm2_to_kkmspc2, kkmspc2_to_wm2),
        (jykms, astrophys.solLum, jykms_to_lsun, lsun_to_jykms),
        (jykms, kkmspc2, jykms_to_kkmspc2, kkmspc2_to_jykms),
        (astrophys.solLum, kkmspc2, lsun_to_kkmspc2, kkmspc2_to_lsun)
    ]
