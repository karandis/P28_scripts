'''
Script to make plots for the paper

Creates the Response correction for the P2 star and the derivation of the response itself
'''

import numpy as np
import pandas as pd
import plotting
from pathlib import Path
import matplotlib.pyplot as pl
import scipy.interpolate as si
from astropy.io import fits
plotting.set_mode('paper')

'''
To use:
307479,V__AG_Psc
'''
unseq=307479
obj = 'V__AG_Psc'

basepath = Path(f'/STER/karansinghd/PhD/Projects/P28_c')

#/STER/karansinghd/PhD/ResponseCorrection/responses_c/20130228_HD46300_453582.txt

# paths= ['/STER/karansinghd/PhD/Projects/P28_c/NSV__6697/corrected/00453599_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits','/STER/karansinghd/PhD/Projects/P28_c/NSV__6687/corrected/00453600_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits','/STER/karansinghd/PhD/Projects/P28_c/NSV__7077/corrected/00453601_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits','/STER/karansinghd/PhD/Projects/P28_c/NAME_ARCTURUS/corrected/00453602_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits','/STER/karansinghd/PhD/Projects/P28_c/___42_Dra/corrected/00453603_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits']

# paths=['/STER/karansinghd/PhD/Projects/P28_c/HD_286340/corrected/00453585_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits','/STER/karansinghd/PhD/Projects/P28_c/NSV__4456/corrected/00453591_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits','/STER/karansinghd/PhD/Projects/P28_c/__gam_Mon/corrected/00453593_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits','/STER/karansinghd/PhD/Projects/P28_c/NSV__4568/corrected/00453595_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits','/STER/karansinghd/PhD/Projects/P28_c/NSV__4456/corrected/00453596_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits','/STER/karansinghd/PhD/Projects/P28_c/__psi_UMa/corrected/00453597_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits']
#
#
# for p in paths:
#     h = fits.open(p)
#     w = h[1].data.field('wave')
#     f= h[1].data.field('flux')
#     pl.plot(w,f)
#     pl.show()
#
# exit()

c_path=basepath.joinpath(f'{obj}/output/00{unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c_TAC.fits')

hdu = fits.open(c_path)
table = hdu[1].data
w = table.field(0)
f = table.field(2)
f_m = table.field(5)
print(table)


cr_path = basepath.joinpath(f'{obj}/corrected/00{unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits')
hdu=fits.open(cr_path)
tab = hdu[1].data
wave = tab.field('wave')
flux = tab.field('flux')
flux_molec = tab.field('flux_molec')



l='/STER/karansinghd/PhD/ResponseCorrection/responses_c/20100930_HD42818_307484.txt'
df = pd.read_csv(l,sep='\t')



left = min(df.wavelength)+10
right = max(df.wavelength) -10
#create knot array
added = np.array([3781,3852,3913,3952])
violet = np.linspace(4000,4340,25)
blue = np.linspace(4350,5515,30)
green = np.linspace(5545,6635,30)
red = np.linspace(6675,right,20)

knots = np.concatenate((added,violet,blue,green,red))


fig, (ax1,ax2,ax3) = pl.subplots(nrows=3,ncols=1,gridspec_kw={'height_ratios':[1,1,1]},figsize=(8,9),sharex=True)
ax1.plot(w,f_m/10000,'k')
ax1.set_ylim([-0.15,6.3])
ax1.set_ylabel('Flux ($ f\,_{obj}^{0,FF}$)')
# ax2.plot(df.wavelength,df.response*100,'r',label=r'$\frac{f\,^{0,FF}_{cal}}{S^0_{CAL}}$')
ax2.plot(df.wavelength,df.spline*100,'k',ls='-',label=r'$\hat{\Re}$')
# [ax2.axvline(x=xline,color='grey',ls='-',alpha=0.5) for xline in knots]
ax2.set_ylim([0.67,2.9])
ax2.set_ylabel('Effective response ($\hat{\Re}$)')
ax3.plot(wave,flux_molec/1e6,'k')
ax3.set_ylabel('Corrected flux ($S^0_{OBJ}$)')
ax3.set_ylim([-0.1,6])
ax3.set_xlim([3893,9020])
# ax2.legend()
ax3.set_xlabel('Wavelength ($\AA$)')
# ax1.plot(wave,flux)
# ax3.plot(wave,r1)
# ax3.plot(wave,p1)
# ax2.plot(w,f)
# ax2.plot(w,p1*f)
# ax2.plot(w,r1,'k')
fig.tight_layout()
# fig.savefig('/Users/karansinghd/PhD/Projects/P28/Correction_HD3379.pdf',dpi=300,format='pdf')

'''
Fig 2, with the P2 star response derivation...
'''
night=20100930
unseq_p2 = 307484
starname = 'HD42818'
sed_path = Path('/STER/karansinghd/PhD/ResponseCorrection/ModelSEDs')

c_path_p2 = Path(f'/STER/karansinghd/PhD/ResponseCorrection/Molecfit/HD42818/output/00{unseq_p2}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c_TAC.fits')

hdu = fits.open(c_path_p2)
table = hdu[1].data
w = table.field(0)
f = table.field(2)
f_m = table.field(5)
# print(table)


ref_wave, ref_flux = np.genfromtxt(fname=f'{sed_path}/{starname}.rgs', unpack=True)

fig1, (ax11,ax22,ax33) = pl.subplots(nrows=3,ncols=1,gridspec_kw={'height_ratios':[1,1,1]},figsize=(8,9),sharex=True)


ax11.plot(w,f_m/1e4,'k')
ax11.set_ylim([-0.3,18])
ax11.set_ylabel('Flux ($ f\,_{cal}^{0,FF}$)')

ax22.plot(ref_wave,ref_flux/1e6,'k')
ax22.set_ylabel(r'Model flux ($S^0_{CAL}$)')
# [ax22.axvline(x=xline,color='grey',ls='-',alpha=0.5) for xline in knots]

# ax22.set_ylim([-0.1,6])


ax33.plot(df.wavelength,df.response*100,'r',label=r'$\frac{f\,^{0,FF}_{cal}}{S^0_{CAL}}$')
ax33.plot(df.wavelength,df.spline*100,'k',ls='--',label=r'$\hat{\Re}$')
[ax33.axvline(x=xline,color='grey',ls='-',alpha=0.5) for xline in knots]
ax33.set_ylim([0.67,2.9])
ax33.set_ylabel('Effective response ($\hat{\Re}$)')
ax33.set_xlabel('Wavelength ($\AA$)')
ax33.set_xlim([3893,9020])

ax33.legend()


# #interpolate the model at the wavelengths of the spectrum
# model_flux = si.interp1d(ref_wave, ref_flux, kind='linear',fill_value='extrapolate')(wave_p2)

fig1.tight_layout()
# fig1.savefig('/Users/karansinghd/PhD/Projects/P28/ResponseDerivation.pdf',dpi=300,format='pdf')

pl.show()

# night=20100930
# for i in range(5):
#     resp_paths=list(Path('/STER/karansinghd/PhD/ResponseCorrection/responses_c/').glob(f'{night}_*.txt'))
#     print(resp_paths[0])
    # continue
