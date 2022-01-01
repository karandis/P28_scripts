import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as pl
import plotting
import scipy.interpolate as si
from astropy.io import fits
import subprocess as sp
import matplotlib as mpl
import shutil
plotting.set_mode('paper')
mpl.rcParams['agg.path.chunksize'] = 10000
destination = Path('/STER/karansinghd/P28_c_cr')

def copy_spectrum(fitsfile):
    filename = fitsfile.name
    outpath = destination.joinpath(filename)
    # print(fitsfile,filename,outpath)
    # exit()
    shutil.copyfile(fitsfile,outpath)
    return

def load_HERMES_data(fitsfile):
    hdu=fits.open(fitsfile)
    head=hdu[0].header
    crpix = head['CRPIX1']-1
    n_points = head['NAXIS1']
    delta_w = head['CDELT1']
    start_w = head['CRVAL1'] - crpix * delta_w
    logw = np.linspace(start_w, start_w+(delta_w * (n_points-1)), n_points)
    wave = np.exp(logw)
    flux=hdu[0].data
    return wave,flux

def load_TAC_file(fitsfile):
    table2 = fits.open(fitsfile)[1].data
    w2 = table2.field(0)
    f3 = table2.field('flux')
    f2 = table2.field(5)
    return w2,f2,f3

def load_corrected_fits(fitsfile):
    hdu = fits.open(fitsfile)
    flag = hdu[0].header['STDNIGHT']
    table = hdu[1].data
    w = table.field('wave')
    f = table.field('flux')
    fm = table.field('flux_molec')
    return w,f,flag,fm

P2_path = Path('/STER/karansinghd/PhD/Projects/P28_p2/')
P2_tac_path = Path('/STER/karansinghd/PhD/ResponseCorrection/Molecfit')

df = pd.read_csv('p2_table.csv',header=0)
# df = pd.read_csv('melchiors_meta_ckcc.txt',header=0,sep='|')
# df.columns = [x.strip() for x in df.columns]
# df.obsid = df.obsid.astype('int32')
# df.pop('Unnamed: 0')
# df.pop('Unnamed: 29')

df = df.drop_duplicates(subset=['stdunseq'])
# print(len(cr_list),len(df))
# exit()
for i,row in df.iterrows():

    night = int(row['stdnight'])
    unseq = int(row['stdunseq'])
    print(night,unseq)
    # if unseq<769700:
    #     continue
    object = row['stdname']

    cr_path = P2_path.joinpath(f'00{unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits')
    tac_path = P2_tac_path.joinpath(f'{object}/output/00{unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c_TAC.fits')
    c_path = Path(f'/STER/mercator/hermes/{night}/reduced/00{unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits')
    outpath = P2_path.joinpath(f'{unseq}_diagnostic.png')

    if cr_path.is_file():
        # print(cr_path)
    #     copy_spectrum(cr_path)
    #     print(cr_path)
    #     continue
        w_cor,f_cor,flag,f_mol = load_corrected_fits(cr_path)
    else:
        # print(cr_path)
        continue
    w_tac,f_tac,f_c = load_TAC_file(tac_path)
    # w_c,f_c = load_HERMES_data(c_path)
    #get ylim for _c and _cr
    idx = (w_tac>=3800) & (w_tac<=3805)
    ulim = np.nanmean(f_c[idx])
    i2 = (w_cor>=3800) & (w_cor<=3805)
    u2 = np.nanmean(f_cor[i2])

    fig,(ax1,ax2) = pl.subplots(nrows=2,ncols=1,gridspec_kw={'height_ratios':[1,1]},figsize=(8,6),sharex=True)
    ax1.plot(w_tac,f_c,'k')
    # ax1.plot(w_tac,f_tac,'r',label='mfit corr')
    ax2.plot(w_cor,f_cor,'k',label='flux')
    ax2.plot(w_cor,f_mol,'r',label='flux_molec')
    ax2.set_xlabel(r'Wavelength ($\AA$)')
    ax1.set_ylabel('Flux')
    ax2.set_ylabel('Flux')
    ax1.set_xlim([3730,9020])
    ax1.set_ylim(bottom=0,top=ulim)
    ax2.set_ylim(bottom=0,top=2*u2)
    ax1.set_title(f'{object}, {night}, {unseq}, STDNIGHT={bool(flag)}')
    fig.tight_layout()
    ax1.grid(ls='--',alpha=0.5)
    ax2.grid(ls='--',alpha=0.5)
    ax2.legend()
    # ax3.grid(ls='--',alpha=0.5)

    fig.savefig(outpath,dpi=300,bbox_inches='tight')
    pl.close()
    # pl.show()
    # exit()
