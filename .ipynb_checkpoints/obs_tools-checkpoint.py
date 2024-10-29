from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import astropy.units as u
import astropy.coordinates as coord
import math
import statistics
import numpy as np



def magnitude_comp(n, positions, flux_count, flux_count_pluto):
    """
    This function computes the magnitude of Pluto 
    based on the flux in counts of a CCD. It first 
    takes the mean value of the k constant in the 
    magnitude equation m = k - 2.5log(Flux), with the values
    of Flux (in counts) and RA DEC of n sources of the CCD 
    and its values of magnitude in the gaia dr3 data.

    n                -----> number of sources used
    positions        -----> vector containing the positions of each source
    it must be in the format '00 00 00.00 00 00 00.00'
    flux_count       -----> vector containing the flux in counts of each source
    taken in the CCD
    flux_count_pluto -----> flux (in counts) of the object we 
    assume that is Pluto
    """
    mags = np.zeros(n)
    k_cte = np.zeros(n)
    for i in range(n):
        radec = positions[i]
        flux = flux_count[i]
        result = Vizier.query_region(coord.SkyCoord(radec, unit=(u.hourangle, u.deg)),
                        radius=Angle(2.0, "arcsec"), 
                        catalog='I/355')
        mags[i] = float(result[0][0]["Gmag"])
        mag = float(result[0][0]["Gmag"])
        kk = mag + 2.5*math.log10(flux)
        k_cte[i]=kk

    k_mean = statistics.mean(k_cte)
    k_std = statistics.stdev(k_cte)
    pluto_mag = k_mean - 2.5*math.log10(flux_count_pluto)

    return k_mean, k_std, pluto_mag
    

        
    
        


    
