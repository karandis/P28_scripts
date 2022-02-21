import numpy as np
import pandas as pd
from pathlib import Path

good_objs =['HD152614','HD36267','HD118098','HD14055','HD87887','HD46300', 'HD184006','HD149212','HD147449','GSC4293-0432','HD214994','HD42818','HD56169','HD185395','HD206826','HD220657','HD84937']

progid  = 28
calibrations = ['CALIBRATION','Wavelength cal','BIAS','Flat','FF','TH','Th']

basedir=Path.cwd()
HermesPath=Path('/STER/mercator/hermes')

#Keep outpath = basedir if you want the output in current working directory,
#otherwise it creates a subfolder
outpath=basedir.joinpath(f'P{progid}')
filelist=outpath.joinpath(f'P{progid}.list')
outpath.mkdir(exist_ok=True)
#################################################
#open the master file with pandas
print('Reading in /STER/mercator/hermes/HermesFullDataOverview.tsv')
data = pd.read_csv("/STER/mercator/hermes/HermesFullDataOverview.tsv",sep='\t',skiprows=2)
print("Done!")
header = ['unseq', 'prog_id', 'obsmode', 'bvcor', 'observer', 'object', 'ra','dec', 'bjd', 'exptime', 'pmtotal', 'date-avg', 'airmass','filename']
data.columns=header

#Search by prog id
data_p2 =data.loc[(data['prog_id']==2) & (~data['object'].str.contains('|'.join(calibrations),na=False))]

#same for p2
nights2 =[data_p2['filename'].values[i].split('/')[4] for i in range(len(data_p2))]
data_p2=data_p2.assign(night=nights2)

data_p2.night=data_p2.night.astype('int32')

data_p2 = data_p2.loc[data_p2.night >=20100927]

data_p2 = data_p2.loc[data_p2.object.isin(good_objs)]

df = pd.read_csv('OverviewProgram28.csv',header=0,sep='|')
df.columns = [x.strip() for x in df.columns]
nights_df =[df['filename'].values[i].split('/')[4] for i in range(len(df))]


# df=df.assign(flag_cr=flag)

df=df.assign(night=nights_df)
df = df.loc[flag]
df.night=df.night.astype('int32')
df = df.loc[df.night >=20100927]

un_n = df.night.unique()
un_p2 =  data_p2.night.unique()
# print(sorted(list(set(un_n) & set(un_p2))))
# print(data_filtered.night.isin(data_p2.night.values).sum())
print(len(list(set(un_n) & set(un_p2))))
print(df.night.isin(data_p2.night.values).sum())
# print(len(data_filtered))
print(len(df.night.unique()))
print(len(df.night))

p2_needed = data_p2.loc[data_p2.night.isin(df.night.unique())]
print(p2_needed.info())
print(len(p2_needed.night.unique()))
print(p2_needed.groupby('object').size())
#NIGHT: 20201116, UNSEQ: 981462, OBJECT: BD+17.4708

exit()
