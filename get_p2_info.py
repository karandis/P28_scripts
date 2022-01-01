import numpy as np
import pandas as pd
import plotting
from pathlib import Path
import matplotlib.pyplot as pl
import scipy.interpolate as si
from astropy.io import fits
import get_response_alhena as gr

#-------------------------------------------------------------------------------

bad_models = ['GSC4293-0432','HD184006','HD147449','HD185395','HD206826','HD84937','HD220657']

# d = pd.read_csv('/STER/karansinghd/P28_c_cr/stdinfo.csv')
# print(d.head())
# exit()
d = pd.read_csv('/STER/karansinghd/PhD/Projects/P28_c/stdinfo.csv')
d.columns = ['unseq','stdunseq','stdname','night','stdnight']
# d = pd.read_csv('/Users/karansinghd/PhD/Projects/P28/p2_table.csv')
d.unseq = d.unseq.astype('int32')
d.stdunseq = d.stdunseq.astype('int32')
d.night = d.night.astype('int32')
d.stdnight = d.stdnight.astype('int32')

print(d.loc[d.stdname.isin(bad_models)])
print(len(d))
d.to_csv('/STER/karansinghd/P28_c_cr/stdinfo.csv',index=False)
exit()
'''
snippet to get objects per night
# d = d.drop_duplicates(subset=['night'])
# print(d.groupby('stdname').size().sum())
'''
n = []
for p in Path('/STER/pierre/hermes/p28/funkyblue').glob('*.png'):
    n.append(int(p.stem.split('_')[0]))

df = d.loc[d.unseq.isin(n)]
# print(d.unseq)
print(df.stdname.unique())
print(df.loc[df.stdname == 'GSC4293-0432'])#374008
# print(df.loc[df.unseq == 374008])#374008
# print(n)
exit()

#-------------------------------------------------------------------------------
resp_path = Path('/STER/karansinghd/PhD/ResponseCorrection/responses_c/')
#Read in P28 meta data
df = pd.read_csv('melchiors_meta_ckcc.txt',header=0,sep='|')
df.columns = [x.strip() for x in df.columns]
df.obsid = df.obsid.astype('int32')
df.pop('Unnamed: 0')
df.pop('Unnamed: 29')

#read the overview file to get the night
dfo = pd.read_csv('/STER/karansinghd/PhD/Projects/P28_c/OverviewProgram28.csv',header=0,sep='|')
dfo.columns = [x.strip() for x in dfo.columns]
nights_dfo =[dfo['filename'].values[i].split('/')[4] for i in range(len(dfo))]

dfo=dfo.assign(night=nights_dfo)
dfo.night=dfo.night.astype('int32')

df_info = pd.DataFrame({'unseq':[],'stdunseq':[],'stdname':[],'night':[],'stdnight':[]})
print(df_info)

nights_needed=[]
for j,row in df.iterrows():
    night = int(dfo.loc[row.obsid == dfo.unseq].night.values[0])
    # print(night)
    # exit()
    nights_needed.append(dfo.loc[row.obsid == dfo.unseq].night.values[0])
    # paths = list(resp_path.glob(f'{night}_*.txt'))
    # if len(paths):
    #     fname = paths[0].name
    #     stdname = fname.split('_')[1]
    #     stdunseq = int(fname.split('_')[2].split('.')[0])
    #     stdnight=int(night)
    # else:
    #     resp = gr.response(night=night, tolerance=90, overwrite = False)
    #     indices = resp.unseq.keys()
    #     for k in indices:
    #         if Path(resp.rpath[k]).is_file():
    #             fname = Path(resp.rpath[k]).name
    #             stdname = fname.split('_')[1]
    #             stdunseq = int(fname.split('_')[2].split('.')[0])
    #             stdnight = int(fname.split('_')[0])
    #             break

    # df_info.loc[j]=[int(row.obsid),stdunseq,stdname,night,stdnight]
    # print(row.obsid,stdunseq,stdname)


# df_info.to_csv('/STER/karansinghd/PhD/Projects/P28_c/p2_table.txt',sep=',',index=False,header=True)
# df_info.to_csv('/Users/karansinghd/PhD/Projects/P28/p2_table.csv',sep=',',index=False,header=True)

df=df.assign(night=nights_needed)

df.to_csv('melchiors_meta.csv',index=False,sep='|')

# night=20100930
# for i in range(5):
#     resp_paths=list(Path('/STER/karansinghd/PhD/ResponseCorrection/responses_c/').glob(f'{night}_*.txt'))
#     print(resp_paths[0])
    # continue
