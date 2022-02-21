from pathlib import Path
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as pl
import plotting
import numpy as np

#Initialise paths
funkyblue_dir = Path('/STER/pierre/hermes/p28/funkyblue')

cf_cr = Path('/STER/karansinghd/P28_crcf')
c_cr = Path('/STER/karansinghd/P28_c_cr')

#load files
stdinfo = pd.read_csv('/STER/karansinghd/PhD/Projects/P28_c/stdinfo.csv',header=0)

stdinfo.columns=['unseq', 'stdunseq', 'stdname', 'night', 'stdnight']
stdinfo=stdinfo.astype({'unseq':int, 'stdunseq':int, 'stdname':str, 'night':int, 'stdnight':int})

melchiors_meta = pd.read_csv('/STER/karansinghd/PhD/Projects/P28_c/melchiors_meta.csv',header=0,sep='|')

'''
tests to be done:
1) make a list of funkyblue
2) compare crcf and ccr to find list where one is there but not the other
3) From 2), isolate STDNIGHT=1/0
'''

fb_list = []

for f in funkyblue_dir.glob('*.png'):
    fb_list.append(f.stem.split('_')[0])

std_list = stdinfo.loc[stdinfo['unseq'].isin(fb_list)]

# print(std_list.loc[std_list.stdname=='HD42818'].stdnight.value_counts())
# print(std_list.loc[std_list.stdnight==20100927])

# print(stdinfo.stdname.value_counts())
#snippet to find SNR from HERMES overview

data = pd.read_csv("/STER/mercator/hermes/HermesFullDataOverview.tsv",sep='\t',skiprows=2)
header = ['unseq', 'prog_id', 'obsmode', 'bvcor', 'observer', 'object', 'ra','dec', 'bjd', 'exptime', 'pmtotal', 'date-avg', 'airmass','filename']
data.columns=header

# print(data.loc[data.unseq.isin(std_list.stdunseq.unique())])

d_fb = data.loc[data.unseq.isin(std_list.stdunseq.unique())]
# d_fb = data.loc[data.unseq.isin(stdinfo.stdunseq.unique())]

d_fb.reset_index(inplace=True, drop=True)
print(d_fb)
# exit()
snv = []

for i,row in d_fb.iterrows():
    if i%10==0:
        print(i)
    ext=fits.open(row['filename'])
    night = row['filename'].split('/')[4]
    head = fits.open(f'/STER/mercator/hermes/{night}/reduced/00{row["unseq"]}_HRF_OBJ_ext_CosmicsRemoved_log_merged_cf.fits')[0].header
    print(row["unseq"],str(head['*SNR*']).strip())
    continue
    exit()
    # usignal = np.median(ext[0].data[2050:2250,51])
    # bsignal = np.median(ext[0].data[2250:2350,41])
    vsignal = np.median(ext[0].data[2350:2450,24])
    # rsignal = np.median(ext[0].data[2450:2650,15])
    # isignal = np.median(ext[0].data[2300:3100,5])
    # print(np.sqrt(vsignal))
    snv.append(np.sqrt(vsignal))

exit()
    # exit()
# d_fb["SNRV"] = snv
# d_fb.to_csv('FunkyBlueData.csv',index=False)
# d_fb = pd.read_csv('FunkyBlueData.csv')

print('done')
# hist = d_fb.hist(column='SNRV')
# print(d_fb.SNRV.min(),d_fb.SNRV.max())
# pl.ylabel('Frequency')
# pl.xlabel('SNR_V')
# pl.savefig('/Users/karansinghd/PhD/Projects/ResponseCorrection/2022/SNRV_all.png')

fig,a = pl.subplots()
d_fb.hist(column='airmass')
pl.ylabel('Frequency')
pl.xlabel('Airmass')
pl.savefig('/Users/karansinghd/PhD/Projects/ResponseCorrection/2022/Airmass.png')
pl.show()
# print(data.head())
