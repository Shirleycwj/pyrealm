import pyreadr
import pyrealm
from pyrealm import pmodel
import numpy as np
import pandas as pd
import pmodel_function
import data_preparation as dp

elev = pd.read_csv("/Users/wenjia/Documents/PhD/ASC_201008_/global_simulation_201203/inputs/ele_map.csv", header=None)
self_elev = np.asarray(elev)
c4_percent = pd.read_csv("/Users/wenjia/OneDrive - Imperial College London/PhD_200803/Research/20190415_C4/c4_percent_map.csv", header=None)
c3_percent = pd.read_csv("/Users/wenjia/OneDrive - Imperial College London/PhD_200803/Research/20190415_C4/c3_percent.csv", header=None)
c3 = np.asarray(c3_percent)
c4 = np.asarray(c4_percent)
co2_monthly = pd.read_csv("/Users/wenjia/Documents/PhD/ASC_201008_/global_simulation_201203/inputs/co2_monthly_01_18.csv", header=None)
self_co2 = np.array(co2_monthly[1], dtype=np.float64)

tmp = dp.rda_read('/Users/wenjia/Documents/PhD/ASC_201008_/global_simulation_201203/inputs/climate_210127/tmp_0119.rda')
tmn = dp.rda_read('/Users/wenjia/Documents/PhD/ASC_201008_/global_simulation_201203/inputs/climate_210127/tmn_0119.rda')
tmx = dp.rda_read('/Users/wenjia/Documents/PhD/ASC_201008_/global_simulation_201203/inputs/climate_210127/tmx_0119.rda')
vap = dp.rda_read('/Users/wenjia/Documents/PhD/ASC_201008_/global_simulation_201203/inputs/climate_210127/vap_0119.rda')
# sw_in = dp.rda_read('/Users/wenjia/Desktop/swin_0101_narm.rda') # na.rm SW_in
# sw_in = dp.rda_read('/Users/wenjia/Documents/PhD/ASC_201008_/global_simulation_201203/inputs/climate_210127/swin_0116_scaled.rda')
sw_in = dp.rda_read('/Users/wenjia/Documents/PhD/ASC_201008_/global_simulation_201203/inputs/climate_210127/swin_0118.rda')
# fapar = dp.rda_read('/Users/wenjia/Documents/PhD/ASC_201008_/global_simulation_201203/inputs/climate_210127/fapar_8216.rda')
fapar = dp.rda_read('/Users/wenjia/Documents/PhD/ASC_201008_/global_simulation_201203/inputs/climate_210127/fapar_0116.rda')

ap = dp.rda_read('/Users/wenjia/Documents/PhD/ASC_201008_/global_simulation_201203/inputs/soilLim_210312/alpha.rda')
sm = dp.rda_read('/Users/wenjia/Documents/PhD/ASC_201008_/global_simulation_201203/inputs/soilLim_210312/soillim.rda')
# soil_lim = dp.rda_read('/Users/wenjia/Documents/PhD/ASC_201008_/global_simulation_201203/inputs/wd_stress_0218_cor.rda')

for year in range(1901, 2017):
    # year = 1982
    gpp_mon = np.zeros(shape=(360, 720, 12))
    for mon in range(0, 12):
        num = (year-1901) * 12 + mon
        # num_2 = (year-1982) * 12 + mon
        # num_fapar = mon  # using 1982 fapar for 1901-1981
        # num_ppfd = (year-1979) * 12 + mon
        # num_sl = (year-2002) * 12 + mon

        # num = 0
        self_tmn = np.asarray(tmn[:, :, num])
        self_tmx = np.asarray(tmx[:, :, num])
        self_tmp = np.asarray(tmp[:, :, num])
        self_vap = np.asarray(vap[:, :, num])
        self_ppfd = np.asarray(sw_in[:, :, num] * 2.04/1000)  # convert radiation to ppfd
        self_fapar = np.asarray(fapar[:, :, num])
        self_ap = np.asarray(ap[:, :, num])
        self_sm = np.asarray(sm[:, :, num])
        # self_soil_lim = np.asarray(soil_lim[:, :, num_sl])

        self_p = pmodel.calc_patm(self_elev)
        self_vpd = pmodel_function.calc_vpd(self_p, self_vap, self_tmx, self_tmn)
        self_vpd[self_vpd < 0] = 0

        model_c3 = pmodel.PModel(tc=self_tmp, patm=self_p, vpd=self_vpd, co2=self_co2[num], soilmstress=1)
        out_lue_c3 = model_c3.unit_iabs.scale_iabs(fapar=self_fapar, ppfd=self_ppfd).lue
        gpp_c3 = model_c3.unit_iabs.scale_iabs(fapar=self_fapar, ppfd=self_ppfd).gpp * c3
        gpp_c3[out_lue_c3 < 0] = 0

        model_c4 = pmodel.PModel(tc=self_tmp, patm=self_p, vpd=self_vpd, co2=self_co2[num], soilmstress=1, c4=True)
        out_lue_c4 = model_c4.unit_iabs.scale_iabs(fapar=self_fapar, ppfd=self_ppfd).lue
        gpp_c4 = model_c4.unit_iabs.scale_iabs(fapar=self_fapar, ppfd=self_ppfd).gpp * c4
        gpp_c4[out_lue_c4 < 0] = 0

        sm_flue = pmodel.calc_soilmstress(self_sm, self_ap)

        gpp = (gpp_c3 + gpp_c4) * sm_flue
        gpp_mon[:, :, mon] = gpp
        np.savetxt('/Users/wenjia/Documents/PhD/ASC_201008_/global_simulation_201203/Ra_simu/global_1901_2016/gpp_mon/'+str(year)+str('%02d' % (mon+1))+'.csv', gpp, delimiter = ',')

    gpp_annual = np.nansum(gpp_mon, axis=2) * 30
    # print(gpp_annual.shape)

    np.savetxt('/Users/wenjia/Documents/PhD/ASC_201008_/global_simulation_201203/Ra_simu/global_1901_2016/gpp_annual/'+str(year)+'.csv', gpp_annual, delimiter=',')




# tmn_mon = tmn[:, :, 1427]
# print(tmn_mon.shape)
# swin_mon = swin[:, :, 0]
# print(swin_mon.shape)
