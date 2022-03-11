import pyreadr
import pyrealm
from pyrealm import pmodel
import numpy as np
import pandas as pd
import data_preparation as dp
from pyrealm.utilities import convert_vp_to_vpd
from pyrealm.param_classes import HygroParams

elev = pd.read_csv("/Volumes/Extreme SSD/Base data/ele_map.csv", header=None)
self_elev = np.asarray(elev)
c4_percent = pd.read_csv("/Volumes/Extreme SSD/Base data/c4_percent_map.csv", header=None)
c3_percent = pd.read_csv("/Volumes/Extreme SSD/Base data/c3_percent.csv", header=None)
c3 = np.asarray(c3_percent)
c4 = np.asarray(c4_percent)
tmp = dp.rda_read('/Volumes/Extreme SSD/Base data/tmp_0119.rda')
tmn = dp.rda_read('/Volumes/Extreme SSD/Base data/tmn_0119.rda')
tmx = dp.rda_read('/Volumes/Extreme SSD/Base data/tmx_0119.rda')
vap = dp.rda_read('/Volumes/Extreme SSD/Base data/vap_0119.rda')
sw_in = dp.rda_read('/Volumes/Extreme SSD/Base data/swin_0118.rda')
fapar = dp.rda_read('/Volumes/Extreme SSD/Base data/fapar_0116.rda')
co2_monthly = pd.read_csv("/Volumes/Extreme SSD/Base data/co2_monthly_01_18.csv", header=None)
self_co2 = np.array(co2_monthly[1], dtype=np.float64)
ap = dp.rda_read('/Volumes/Extreme SSD/Base data/alpha.rda')  # AET/PET
sm = dp.rda_read('/Volumes/Extreme SSD/Base data/soillim.rda')  # soil moisture

for year in range(1982, 2017):
    gpp_mon = np.zeros(shape=(360, 720, 12))
    for mon in range(0, 12):
        num = (year-1901) * 12 + mon

        self_tmn = np.asarray(tmn[:, :, num])
        self_tmn[self_tmn <= -25] = float('nan')
        self_tmx = np.asarray(tmx[:, :, num])
        self_tmx[self_tmx <= -25] = float('nan')
        self_tmp = np.asarray(tmp[:, :, num])
        self_tmp[self_tmp <= -25] = float('nan')

        self_ppfd = np.asarray(sw_in[:, :, num] * 2.04/1000)  # convert radiation to ppfd
        self_fapar = np.asarray(fapar[:, :, num])
        self_ap = np.asarray(ap[:, :, num])
        self_sm = np.asarray(sm[:, :, num])

        self_vap = np.asarray(vap[:, :, num])
        self_p = pmodel.calc_patm(self_elev)
        self_vap_elev = self_vap * (self_p/pmodel.calc_patm(0))  # calculate vap at certain altitude in hPa

        hygro_par = HygroParams(magnus_option='Allen1998')  # calculating VPD using parameters from Allen1998
        vpd_min = convert_vp_to_vpd(vp=0.1*self_vap_elev, ta=self_tmn, hygro_params=hygro_par)
        vpd_max = convert_vp_to_vpd(vp=0.1*self_vap_elev, ta=self_tmx, hygro_params=hygro_par)
        vpd_rng = np.stack([vpd_min, vpd_max])
        self_vpd = vpd_rng.mean(axis=0) * 1000  # convert to Pa
        self_vpd[self_vpd < 0] = 0

        env = pmodel.PModelEnvironment(tc=self_tmp, co2=self_co2[num], patm=self_p, vpd=self_vpd) # calculating environment
        sm_flue = pmodel.calc_soilmstress(soilm=self_sm, meanalpha=self_ap)  # estimating soil moisture stress factor

        # Estimate GPP for C3
        model_c3 = pmodel.PModel(env)
        model_c3.estimate_productivity(fapar=self_fapar, ppfd=self_ppfd)
        gpp_c3_weighted = model_c3.gpp * c3
        gpp_c3_weighted[self_tmp < 0] = 0  # force GPP to 0 when temperature below freezing
        gpp_c3_weighted[gpp_c3_weighted < 0] = 0  # force negative GPP to 0
        # Estimate GPP for C4
        model_c4 = pmodel.PModel(env, c4=True)
        model_c4.estimate_productivity(fapar=self_fapar, ppfd=self_ppfd)
        gpp_c4_weighted = model_c4.gpp * c4
        gpp_c4_weighted[self_tmp < 0] = 0
        gpp_c4_weighted[gpp_c4_weighted < 0] = 0
        # monthly daily mean GPP (gC/m2 day)
        gpp_mon[:, :, mon] = (gpp_c3_weighted + gpp_c4_weighted) * sm_flue

    # annual total GPP (gC/m2 yr)
    gpp_annual = np.nansum(gpp_mon, axis=2) * 30
    np.savetxt('/Users/wenjiacai/Desktop/annual_GPP/'+str(year)+'.csv', gpp_annual, delimiter=',')
