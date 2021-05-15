import numpy as np
import pandas as pd
import math


es0 = 0.611

def calc_vpd(p,vap,Tmax,Tmin):
    vap_elev = vap * (p/101325)
    es1 = es0 * np.exp((17.27 * Tmax)/(Tmax+273.3))
    vpd1_hpa = (es1 - 0.1 * vap_elev) * 10  # kPa --> hPa
    es2 = es0 * np.exp((17.27 * Tmin)/(Tmin+273.3))
    vpd2_hpa = (es2 - 0.1 * vap_elev) * 10  # kPa --> hPa
    vpd_hpa = (vpd1_hpa + vpd2_hpa)/2  # calculating vpd by average the two
    vpd_pa = vpd_hpa * 100
    return vpd_pa
