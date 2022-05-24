# Functions to get T and EM info out of the GOES/XRS data
# Based off of https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_tem_calc.pro
# and specifically the CHIANTI versions via
# https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_get_chianti_temp.pro
# and
# https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_get_chianti_em.pro
# 
# Response calculated in sswidl:
#   https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/buildresponse/goes_chianti_response.pro
# the response fits is always called "goes_chianti_response_latest.fits" now
# Version have here is from 21-05-2022 and v10 CHIANTI:
#   https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_chianti_response_latest.fits
# Older V9 CHIANTI one also in repo "goes_chianti_resp_20200812.fits", but have wrong factor in GOES13/15 short
#  
# 
# Values returned << 1% out from sswidl versions, probably due to slight difference in the spline interpolation
# 
# 24-10-2021    IGH    Created
# 25-10-2021    IGH    Updated to return newer CHIANTI v9 responses, and option of abundance and sat
# 23-05-2022    IGH    Updated with new response (CHIANTI v10 and G13/15 short fix)
#                      Default loads new response, but old_ver=True gives 20200812 version
# -----------------------------
import numpy as np
from scipy import interpolate
from astropy.io import fits
# -----------------------------
# -----------------------------
def get_resps(sat=15,cor_not_pho=True,old_ver=False):
#   Returns the the GOES/XRS temperature response functions for long and short channels
# 
#   Response calculated in sswidl:
#       https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/buildresponse/goes_chianti_response.pro
#   the response fits is always called "goes_chianti_response_latest.fits" now
#   Version have here is from 21-05-2022 and v10 CHIANTI:
#       https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_chianti_response_latest.fits
# 
#   Input:
#       sat - Which GOES satellite to use? (default 15)
#       cor_not_pho: - Use coronal not photospheric abundances (default True)
#       old_ver - Use previous responses of 20200812 + wrong G13/15 short: (default False) - just for testing
#   Output:
#       resps[:,0] - Temp resp for 1-8\AA in units of 10^{-55} W m^{-2} cm^{3}
#       resps[:,1] - Temp resp for 0.5-4\AA in units of 10^{-55} W m^{-2} cm^{3}
#       resptmk    - Temp binning of TR ratio in MK

    if old_ver:
        rfile='goes_chianti_resp_20200812.fits'
    else:
        rfile='goes_chianti_response_latest.fits'
    hdulist = fits.open(rfile)
    respdat=hdulist[1].data
    hdulist.close()

    resptmk=np.array(respdat["TEMP_MK"][sat-1])
    
    if cor_not_pho:
        abdun="COR"
    else:
        abdun="PHO"  
    resps=np.empty((101,2))
    resps[:,0]=respdat["FLONG_"+abdun][sat-1]
    resps[:,1]=respdat["FSHORT_"+abdun][sat-1]
   
    return resps, resptmk
# -----------------------------
# -----------------------------
def get_tem(fl,fs,sat=15,cor_not_pho=True,old_ver=False):
    
#   Returns the T and EM for ratio of fluxes in short/long of GOES/XRS channels
# 
#   This comes from https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_tem_calc.pro
#   Response calculated in sswidl:
#       https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/buildresponse/goes_chianti_response.pro
#   the response fits is always called "goes_chianti_response_latest.fits" now
#   Version have here is from 21-05-2022 and v10 CHIANTI:
#       https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_chianti_response_latest.fits
#
#   Input:
#       fl - GOES/XRS flux in long channel, 1-8\AA no scaling, Wm^{-2}
#       fs - GOES/XRS flux in short channel,0.5-4\AA scaling, Wm^{-2}
#            Single number, or array
#       sat - Which GOES satellite to use? (default 15)
#       cor_not_pho: - Use coronal not photospheric abundances (default True)
#   Output:
#       TMK - Temperature in MK
#       EM  - Emission Measure in cm^{-3}
# 
#   NOTE -  No longer need scaling for GOES16/17 as was correct but GOES13/15 short resp wrong (now fixed)
# 

#   Get the TR to work out the EM
    resps, resptmk=get_resps(sat,cor_not_pho,old_ver)

#   Then calculate the resps ratio
    resprat=resps[:,1]/resps[:,0]
    
#   Use scipy cubic spline interpolation
    rat_func=interpolate.interp1d(resprat, resptmk,kind='cubic')
    
#   Ratio can't be outside of values of TR ratio, so
    gfs=np.array(fs)
    gfl=np.array(fl)
    grat=np.array(gfs/gfl)
    grat[grat < resprat.min()]=resprat.min()
    grat[grat > resprat.max()]=resprat.max()
#   Work out the temperature 
    tmk=np.array(rat_func(grat))
    
#   Need the actual TR to work out the EM
    resps, resptmk=get_resps(sat,cor_not_pho)
#   Use scipy cubic spline interpolation
    tr18_func=interpolate.interp1d(resptmk,resps[:,0],kind='cubic')
#   Can't use TMK values at/outside the range, so
#   (might this causes issues - better way of doing this??)
    tmk[tmk < resptmk.min()]=1.0001
    tmk[tmk >resptmk.max()]=resptmk.max()
#   Work out the Emission Measure
    em=np.array(1e55*gfl/tr18_func(tmk))
    em[em < 0]=0.
    em[tmk <= 1.0001]=0.
    
    return tmk, em 
# -----------------------------
# -----------------------------
def get_resprat(sat=15,cor_not_pho=True,old_ver=False):

#   Response calculated in sswidl:
#       https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/buildresponse/goes_chianti_response.pro
#   the response fits is always called "goes_chianti_response_latest.fits" now
#   Version have here is from 21-05-2022 and v10 CHIANTI:
#       https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_chianti_response_latest.fits
# 
#   Input :
#       sat - Which GOES satellite to use? (default 15)
#       cor_not_pho: - Use coronal not photospheric abundances (default True)
#       old_ver - Use previous responses of 20200812 + wrong G13/15: (default False) - just for testing
#   Output:
#       resprat - Ratio of TR_(0.5-4) / TR_(1-8)
#       resptmk - Temp binning of TR ratio in MK
# 
    
    if old_ver:
        rfile='goes_chianti_resp_20200812.fits'
    else:
        rfile='goes_chianti_response_latest.fits'
    hdulist = fits.open(rfile)
    respdat=hdulist[1].data
    hdulist.close()
    
    resptmk=np.array(respdat["TEMP_MK"][sat-1])
    
    if cor_not_pho:
        abdun="COR"
    else:
        abdun="PHO"   
    resprat=np.empty((101))
    resprat=respdat["FSHORT_"+abdun][sat-1]/respdat["FLONG_"+abdun][sat-1]
    
    return resprat, resptmk
# -------------------------------
