# Functions to get T and EM info out of the GOES/XRS data
# Based off of https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_tem_calc.pro
# and specifically the CHIANTI versions via
# https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_get_chianti_temp.pro
# and
# https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_get_chianti_em.pro
# 
# At the moment only doing GOES15 and coronal abundances
# 
# Values returned << 1% out from sswidl versions, probably due to slight difference in the spline interpolation
# 
# 24-10-2021    IGH    Created
# -----------------------------
import numpy as np
from scipy import interpolate

# -----------------------------
def get_resprat():

#   Returns the ratio of the GOES/XRS temperature response in the short to long channels
#   i.e. TR_(0.5-4) / TR_(1-8)
# 
#   This comes from https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_get_chianti_temp.pro
#   At the moment just for GOES15 and coronal -> add options to SATnum and coronal/photospheric in future?
# 
#   Output:
#       resprat - Ratio of TR_(0.5-4) / TR_(1-8)
#       resptmk - Temp binning of TR ratio in MK
# 
    
    resprat=np.array([2.43e-06,3.63e-06,5.31e-06,7.69e-06,1.11e-05,1.59e-05,2.29e-05,3.30e-05,\
       4.74e-05,6.80e-05,9.68e-05,1.36e-04,1.88e-04,2.53e-04,3.33e-04,4.29e-04,5.43e-04,\
       6.78e-04,8.39e-04,1.03e-03,1.26e-03,1.53e-03,1.84e-03,2.21e-03,2.63e-03,3.10e-03,\
       3.63e-03,4.22e-03,4.89e-03,5.62e-03,6.45e-03,7.36e-03,8.38e-03,9.51e-03,1.08e-02,\
       1.22e-02,1.37e-02,1.54e-02,1.74e-02,1.95e-02,2.19e-02,2.45e-02,2.73e-02,3.05e-02,\
       3.40e-02,3.78e-02,4.19e-02,4.65e-02,5.14e-02,5.69e-02,6.30e-02,6.97e-02,7.72e-02,\
       8.57e-02,9.53e-02,1.06e-01,1.18e-01,1.32e-01,1.47e-01,1.63e-01,1.80e-01,1.98e-01,\
       2.17e-01,2.36e-01,2.56e-01,2.76e-01,2.96e-01,3.15e-01,3.35e-01,3.54e-01,3.72e-01,\
       3.90e-01,4.08e-01,4.25e-01,4.42e-01,4.58e-01,4.73e-01,4.88e-01,5.02e-01,5.16e-01,\
       5.29e-01,5.41e-01,5.53e-01,5.64e-01,5.75e-01,5.86e-01,5.96e-01,6.05e-01,6.14e-01,\
       6.23e-01,6.31e-01,6.39e-01,6.46e-01,6.53e-01,6.60e-01,6.66e-01,6.72e-01,6.78e-01,\
       6.84e-01,6.89e-01,6.94e-01])
    resptmk=10**(np.arange(101)*0.02)
    
    return resprat, resptmk
# -------------------------------

# -----------------------------
def get_resps():
#   Returns the the GOES/XRS temperature response functions for long and short channels
# 
#   This comes from https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_get_chianti_em.pro
#   At the moment just for GOES15 and coronal -> add options to SATnum and coronal/photospheric in future?
# 
#   Output:
#       resps[:,0] - Temp resp for 1-8\AA in units of 10^{-55} W m^{-2} cm^{3}
#       resps[:,1] - Temp resp for 0.5-4\AA in units of 10^{-55} W m^{-2} cm^{3}
#       resptmk    - Temp binning of TR ratio in MK
    resps=np.empty((101,2))

    resps[:,0]=np.array([6.27e-05,1.15e-04,2.09e-04,3.68e-04,6.29e-04,1.05e-03,1.69e-03,2.65e-03,\
        4.08e-03,6.16e-03,9.17e-03,1.35e-02,1.95e-02,2.81e-02,3.99e-02,5.61e-02,7.78e-02,\
        1.06e-01,1.43e-01,1.91e-01,2.50e-01,3.23e-01,4.11e-01,5.18e-01,6.46e-01,7.98e-01,\
        9.76e-01,1.18e+00,1.42e+00,1.70e+00,2.02e+00,2.38e+00,2.79e+00,3.24e+00,3.75e+00,\
        4.30e+00,4.92e+00,5.58e+00,6.30e+00,7.07e+00,7.90e+00,8.77e+00,9.68e+00,1.06e+01,\
        1.16e+01,1.26e+01,1.37e+01,1.47e+01,1.58e+01,1.69e+01,1.80e+01,1.91e+01,2.02e+01,\
        2.13e+01,2.23e+01,2.34e+01,2.43e+01,2.52e+01,2.61e+01,2.69e+01,2.77e+01,2.85e+01,\
        2.94e+01,3.02e+01,3.11e+01,3.20e+01,3.30e+01,3.40e+01,3.51e+01,3.62e+01,3.74e+01,\
        3.87e+01,4.00e+01,4.14e+01,4.28e+01,4.43e+01,4.58e+01,4.73e+01,4.89e+01,5.05e+01,\
        5.20e+01,5.36e+01,5.52e+01,5.68e+01,5.84e+01,5.99e+01,6.14e+01,6.28e+01,6.42e+01,\
        6.56e+01,6.69e+01,6.81e+01,6.92e+01,7.03e+01,7.12e+01,7.21e+01,7.30e+01,7.37e+01,\
        7.44e+01,7.50e+01,7.55e+01])
    
    resprat, resptmk = get_resprat()
#   Calculate short channel response from long channel and ratio of short/long, so short = rat * long
    resps[:,1]=resprat * resps[:,0]
   
    return resps, resptmk
# -----------------------------
# -----------------------------
def get_tem(fl,fs):
    
#   Returns the T and EM for ratio of fluxes in short/long of GOES/XRS channels
# 
#   This comes from https://hesperia.gsfc.nasa.gov/ssw/gen/idl/synoptic/goes/goes_tem_calc.pro
#   At the moment just for GOES15 and coronal -> add options to SATnum and coronal/photospheric in future?
#
#   Input:
#       fl - GOES/XRS flux in long channel, 1-8\AA no scaling, Wm^{-2}
#       fs - GOES/XRS flux in short channel,0.5-4\AA scaling, Wm^{-2}
#            Single number, or array
#   Output:
#       TMK - Temperature in MK
#       EM  - Emission Measure in cm^{-3}
# 

# First get the resps ratio
    resprat, resptmk = get_resprat()
    
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
    resps, resptmk=get_resps()
#   Use scipy cubic spline interpolation
    tr18_func=interpolate.interp1d(resptmk,resps[:,0],kind='cubic')
#   Can't use TMK values at/outside the range, so
#   (might this causes issues - better way of doing this??)
    tmk[tmk < resptmk.min()]=1.0001
    tmk[tmk >resptmk.max()]=resptmk.max()
#   Work out the Emission Measure
    em=np.array(1e55*gfl/tr18_func(tmk))
    em[em<0]=0.
    
    return tmk, em 
# -----------------------------