#Institution: Brazilian National Observatory
#Subject: Observational Astronomy
#Author: Mariane Dias de Souza Gomes
#Date: September/2024
#Scripts 

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
                        catalog='I/355')              #GaiaDR3
        mags[i] = float(result[0][0]["Gmag"])
        mag = float(result[0][0]["Gmag"])
        kk = mag + 2.5*math.log10(flux)
        k_cte[i]=kk
    
    #mags_std = statistics.stdev(mags)
    k_mean = statistics.mean(k_cte)
    k_std = statistics.stdev(k_cte)
    pluto_mag = k_mean - 2.5*math.log10(flux_count_pluto)

    return k_mean, k_std, pluto_mag#, mags_std
    


def astrometric_correction(n, positions, year_sci, month_sci):
    """
    This function makes a correction in the ra and dec
    data of the catalog gaia dr3 based on the epoch of
    the CCD image used. It uses the equations:
    α_C = α_0 + μ_α*Δt ; δ_C = δ_0 + μ_δ*Δt, where α_C
    is the correction, α_0 is the right ascention ICRS,
    μ_α is the proper motion in right ascention multiplied
    by a cos(δ) factor, Δt = (year_ccdimage - year_catalog), 
    the year_ccdimage can be obtained in its header. The same
    is valid for the declination correction δ_C, except that
    in its case it's not necessary to multiply the cos(δ) factor.
    The year_ccdimage is in decimal, because it accounts for 
    the month of the observation.
    

    n                -----> number of sources used
    positions        -----> vector containing the positions of each source (CCD)
    it must be in the format '00 00 00.00 00 00 00.00'
    year_sci         -----> the year of the observation  (start of frame exposure)
    month_sci         -----> the year of the observation (start of frame exposure)
    """
    
    year_cat = 2016           #gaia epoc     
    month_sci_decim = round(30*month_sci/365,2)
    deltat = (year_sci+month_sci_decim - year_cat)
    
    ra_corr = np.zeros(n)
    de_corr = np.zeros(n)
    delta_ra = np.zeros(n)
    delta_de = np.zeros(n)   
    
    for i in range(n):
        radec = positions[i]
        result = Vizier.query_region(coord.SkyCoord(radec, unit=(u.hourangle, u.deg)),
                            radius=Angle(2.0, "arcsec"), 
                            catalog='I/355')
        ra = float(result[0][0]["RA_ICRS"])
        de = float(result[0][0]["DE_ICRS"])
        pmra_mas = float(result[0][0]["pmRA"])
        pmra_deg = (pmra_mas*1e-3)/3600
        pmde_mas = float(result[0][0]["pmDE"])
        pmde_deg = (pmde_mas*1e-3)/3600

        ra_corr[i] = ra + (pmra_deg/np.cos(de))*deltat
        de_corr[i] = de + pmde_deg*deltat
        
        
    return ra_corr, de_corr

def SNR(gain, rad, flux_star, flux_sky):
    """
    This function computes the signal-to-noise ratio based on aperture photometry.
    gain      ----> gain of the ccd
    rad       ----> aperture radius
    flux_star ----> the star flux in counts of the ccd
    flux_sky  ----> the sky flux in counts 
    """
    SNR = (flux_star*gain)/np.sqrt(flux_star*gain + flux_sky*gain*np.pi*rad**2)
    return SNR