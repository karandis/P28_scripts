#---------------------------------------------------------------------------
#Modified to input files from Pierre's ccf
#---------------------------------------------------------------------------
import get_response_alhena as gr
import pandas as pd
from astropy.io import fits
from scipy import interpolate as si
from pathlib import Path
import numpy as np
import matplotlib.pyplot as pl
import matplotlib as mpl
import os
import subprocess as sp
import logging
import traceback
from scipy.interpolate import LSQUnivariateSpline as us
from astropy.time import Time
import shutil
#---------------------------------------------------------------------------
HermesPath = Path('/STER/mercator/hermes/')
LocalPath = Path('/STER/karansinghd/PhD/Projects/P28_c/')
Mfit_dir = Path('/STER/karansinghd/PhD/ResponseCorrection/Molecfit/')
print(Mfit_dir)
#----------------------------------------------------------------
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s: [%(name)s:%(levelname)s] %(message)s',datefmt='%d-%m %H:%M')
logger = logging.getLogger("CORRECT SPECTRUM")
#----------------------------------------------------------------
bad_models = ['GSC4293-0432','HD184006','HD147449','HD185395','HD206826','HD84937','HD220657']

order_to_be = {'HD36267': 1,
           'HD152614':2,
           'HD149212': 3,
           'HD14055':4,
           'HD87887':5,
           'HD42818': 6,
           'HD214994': 7,
           'HD46300': 8,
           'HD118098': 9,
           'GSC4293-0432': 10,
           'HD56169': 11,
           'HD184006': 12,
           'HD147449': 13,
           'HD185395': 14,
           'HD206826': 15,
           'HD84937': 16,
           'HD220657': 17}

df_stdinfo = pd.DataFrame({'unseq':[],'stdunseq':[],'stdname':[],'night':[],'stdnight':[]})
##########################################################################
class correct_spectrum:
    def __init__(self,night,unseq,object):
        overwrite = True
        self.night = night
        self.unseq = unseq
        self.object = object
        self.stdname = ''
        self.stdunseq = 0
        self.stdnight = 0
        logger.info(f'Processing night: {self.night}, UNSEQ: {self.unseq}, object: {self.object}')
        LocalPath.joinpath(f'{self.object}/corrected').mkdir(exist_ok=True,parents=True)
        self.outpath = LocalPath.joinpath(f'{self.object}/corrected/00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits')
        if self.outpath.is_file() and not overwrite:
            logger.info(f'Spectrum already corrected! {self.outpath}')
            return
        self.filename_var = f'/STER/mercator/hermes/{self.night}/reduced/00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_mergedVar_c.fits'
        self.mfit_object_path = LocalPath.joinpath(self.object)
        self.par_path=self.mfit_object_path/Path(f'{self.object}_{self.unseq}.par')
        # self.filename=Path(f'{self.mfit_object_path}/00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_cf_TAC.fits')
        self.filename=Path(f'{self.mfit_object_path}/output/00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c_TAC.fits')
        # self.header = fits.getheader(f'/STER/mercator/hermes/{self.night}/reduced/00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_cf.fits')
        if not self.filename.is_file():
            logger.info('Molecfit will be run as TAC file does not exist')
            self.prepare_molecfit_files()
            self.modify_par_file()
            try:
                self.rewrite_fits()
            except Exception as e:
                logger.critical('Failed while rewriting FITS file')
                logger.critical(f'Exception:{e}')
                logger.critical(traceback.format_exc())
            self.run_molecfit()
        self.wave, self.flux, self.flux_err, self.flux_no_TAC, self.flux_no_TAC_err,self.header = self.load_FITS_table_data()

        logger.info('Successfully loaded Molecfit corrected spectrum')
        self.response=self.load_response()
        # return

        #For _c
        # idx = (self.wave > 4000)
        # self.wave = self.wave[idx]
        # self.flux = self.flux[idx]
        # self.flux_err = self.flux_err[idx]
        # self.response = self.response[idx]
        self.rc_flux = self.flux/self.response
        self.rc_flux_no_TAC = self.flux_no_TAC/self.response

        self.rc_err = self.flux_err/self.response
        self.rc_err_no_TAC = self.flux_no_TAC_err/self.response
        # f,(a,b,c) = pl.subplots(nrows=3,sharex=True,gridspec_kw={'height_ratios':[1,1,1]})
        # a.plot(self.wave,self.flux)
        # b.plot(self.wave,self.response)
        # c.plot(self.wave,self.rc_flux)
        # pl.show()
        # exit()
        # self.create_specific_plot(self.wave,self.rc_flux,linewidth=2,xlabel='Wavelengh ($\AA$)',ylabel='Response Corrected Flux',title=f'{self.object}_{self.unseq}')
        logger.info(f'STDNAME: {self.stdname}, STDNIGHT: {self.stdnight}, STDUNUSEQ: {self.stdunseq}')
        self.write_fits()



    def load_response(self):
        resp_paths=list(Path('/STER/karansinghd/PhD/ResponseCorrection/responses_c/').glob(f'{self.night}_*.txt'))
        if len(resp_paths):
            self.STD_NIGHT=1
            logger.info(f'{len(resp_paths)} responses available for {self.night}')

            #
            # order = [x.stem.split('_')[1] in bad_models for x in resp_paths]

            resp_paths_sorted = sorted(resp_paths, key=lambda x: order_to_be[x.stem.split('_')[1]])

            # res = sorted(zipped, key = lambda x: x[1])
            # order = [x.stem.split('_')[1] in bad_models for x in resp_paths]
            # resp_paths_sorted = [x for _,x in sorted(zip(order,resp_paths))]

            for p in resp_paths_sorted:
            # for p in resp_paths:
                # print(p.stem.split('_')[1])
                # exit()
                logger.info(f'Processing response: {p}')
                self.stdname=p.stem.split('_')[1]
                self.stdunseq=int(p.stem.split('_')[2])
                self.stdnight=int(p.stem.split('_')[0])
                resp,poly = self.load_and_interpolate(p)
                self.create_plots(resp,poly)
                return poly
                # use = input('Do you want to correct with this response? ([y]/n)\n')
                # if use.lower() == 'n':
                    # continue
                # else:
                    # return poly
        else:
            self.STD_NIGHT = 0
            logger.info('get_response.py has to be run!')
            # t = np.abs(int(input('How many days (int) tolerance should be considered?\n')))
            # if t == 0:
            #     t = None
            resp = gr.response(night=self.night, tolerance=60, overwrite = False)
            logger.info(f'{len(resp.unseq)} responses available to use')

            # if resp:
            #     logger.info(f'{len(resp.unseq)} responses available to use')
            # else:
            #     logger.critical(f'No responses available within {t} days of {self.night}, sorry!')
            #     exit()
            indices = resp.unseq.keys()
            resp_paths = []
            for j in indices:
                try:
                    resp_paths.append(Path(resp.rpath[j]))
                except:
                    continue

            resp_paths_sorted = sorted(resp_paths, key=lambda x: order_to_be[x.stem.split('_')[1]])

            # order = [x.stem.split('_')[1] in bad_models for x in resp_paths]
            # resp_paths_sorted = [x for _,x in sorted(zip(order,resp_paths))]

            for p in resp_paths_sorted:
            # for j in indices:
                if p.stem.split('_')[1] in bad_models:
                    continue
                try:
                    resp,poly = self.load_and_interpolate(p)
                    self.stdname=p.stem.split('_')[1]
                    self.stdunseq=int(p.stem.split('_')[2])
                    self.stdnight=int(p.stem.split('_')[0])
                    self.create_plots(resp,poly)
                    return poly
                except:
                    continue
                # use = input('Do you want to correct with this response? ([y]/n)\n')
                # if use.lower() == 'n':
                #     continue
                # else:
                #     return poly

    def load_and_interpolate(self,p):
        data = pd.read_csv(p,sep='\t',header=0,encoding='utf-8',engine='python')

        idx1=(self.wave>=min(data['wavelength'])) & (self.wave<=max(data['wavelength']))

        self.wave = self.wave[idx1]
        self.flux = self.flux[idx1]
        self.flux_no_TAC = self.flux_no_TAC[idx1]
        self.flux_err = self.flux_err[idx1]
        self.flux_no_TAC_err = self.flux_no_TAC_err[idx1]

        r1 = si.interp1d(data['wavelength'],data['response'],kind='linear')(self.wave)
        try:
            p1 = si.interp1d(data['wavelength'],data['spline'],kind='linear')(self.wave)
        except:
            p1 = self.fit_spline(r1)
        return r1,p1

    def get_atm_ext_cor(self, ext_wl, aext_coeff, head, wl):
        m = (min(wl) < ext_wl) & (ext_wl < max(wl))
        ref_aext = si.interp1d(ext_wl[m], aext_coeff[m], kind='quadratic', bounds_error=False, fill_value='extrapolate')
        z = np.deg2rad(90. - float(head['TELALT'])) # Zenithal angle in radians
        logger.debug(f'Zenithal distance: {z:2.5f} rad')
        atm_cor = 10**(0.4*ref_aext(wl)/np.cos(z))
        return atm_cor

    def clean_nan(self,x0, y0, z0, x1):

        m = (y0 != y0) | (z0 != z0)

        x = x0[~m]
        y = y0[~m]
        z = z0[~m]
        x2= x1[~m]

        try:
            y = si.interp1d(x, y, kind='linear', bounds_error=True)(x0)
            z = si.interp1d(x, z, kind='linear', bounds_error=True)(x0)
            x = x0
            x2=x1
        except ValueError:
            pass

        return x, y, z, x2

    def clean_zeros(self,x0, y0, z0,x1):

        m = (y0 == 0) | (z0 == 0)

        x = x0[~m]
        x2 = x1[~m]
        y = y0[~m]
        z = z0[~m]

        y = si.interp1d(x, y, kind='linear', fill_value='extrapolate')(x0)
        z = si.interp1d(x, z, kind='linear', fill_value='extrapolate')(x0)

        return x0, y, z, x1

    def load_FITS_table_data(self):
        '''
        Function to load table data from the FITS file output by Molecfit.
        '''
        hdul = fits.open(self.filename)
        header = hdul[0].header
        table_data = hdul[1].data
        wave = table_data.field(0)
        flux_err = table_data.field(6)
        flux_TAC = table_data.field(5)
        flux_no_TAC = table_data.field(2)
        flux_no_TAC_err = table_data.field(3)
        return wave, flux_TAC, flux_err, flux_no_TAC, flux_no_TAC_err,header

    def load_HERMES_data(self):
        ''' Function to load a HERMES spectrum.
        Input:
            unseq : Unique sequence that identifies the object.
            @type unseq: int
            night : yyyymmdd
            @type night: int
        Returns:
            wavelength : numpy array of wavelengths(lambda or ln(lambda)) in Angstrom
            flux : numpy array with flux (photons?)
            var : variance of the flux
            head : header of the fits file
        '''
        #The input file names, modified to take them from Pierre's directories

        try:
            flux, head = fits.getdata(f'/STER/mercator/hermes/{self.night}/reduced/00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits',header=True)
        except:
            logger.critical("check filename Pierre")

        try:
            var = fits.getdata(f'/STER/mercator/hermes/{self.night}/reduced/00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_mergedVar_c.fits')
        except:
            try:
                var = fits.getdata(f'/STER/mercator/hermes/{self.night}/reduced/{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_mergedVar_c.fits')
                self.filename_var = f'/STER/mercator/hermes/{self.night}/reduced/{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_mergedVar_c.fits'
            except:
                logger.critical(f'Var filename is wrong! Please check: rewrite_fits()')
                return

        crpix = head['CRPIX1']-1
        n_points = head['NAXIS1']
        delta_w = head['CDELT1']
        start_w = head['CRVAL1'] - crpix * delta_w

        offset = self.correct_bar_vel_shift(head)
        start_w_corr = start_w+ offset
        log_wl_original = np.linspace(start_w, start_w+(delta_w * (n_points-1)), n_points)
        log_wl = np.linspace(start_w_corr, start_w_corr+(delta_w * (n_points-1)), n_points)
        wave = np.exp(log_wl)
        wave_original = np.exp(log_wl_original)

        #Check for and clean 1) nans and 2) zeros
        wave, flux, var,wave_original = self.clean_nan(wave,flux,var,wave_original)
        wave, flux, var,wave_original = self.clean_zeros(wave,flux,var,wave_original)
        logger.debug('Cleaned NaNs and zeros')

        return wave, flux, var, head, wave_original

    def correct_bar_vel_shift(self,Head):
        '''Function to undo barycentric correction.
        Input:
            Header of fits (to get BVCOR): Km/s
        Returns:
            offset
        '''
        if 'BVCOR' in Head:
            bar_v = Head['BVCOR'] #in Km/s
        elif 'VHELIO' in Head:
            bar_v = Head['VHELIO']
        else:
            logger.warning('BVCOR or VHELIO not found, set to 0')
            bar_v=0
        #Doppler shift corrected
        logger.info(f'Barycentric velocity from header: {bar_v} km/s')
        offset = np.log(1-(bar_v)/3.0e5)
        return offset

    def prepare_molecfit_files(self):
        #create directories
        if not os.path.isdir(self.mfit_object_path):
            os.mkdir(self.mfit_object_path)
        if not os.path.isdir(self.mfit_object_path/Path('output')):
            os.mkdir(self.mfit_object_path/Path('output'))
        #copy necessary parameter files
        # Command1 = sp.run(f'cp {self.Mfit_dir}/template.par {self.par_path}', shell=True)
        # print(f'Command {Command1.args} ran with returncode {Command1.returncode}')string = f'I am {num:{".2f" if ppl else ""}}'
        Command2 = sp.run(f'cp {Mfit_dir}/include.dat {self.mfit_object_path}/include.dat', shell=True)
        logger.debug(f'Command {Command2.args} ran {"successfully" if Command2.returncode ==0 else "unsuccessfully"}')
        # logger.debug(f'Command {Command2.args} ran with returncode {Command2.returncode}')
        Command3 = sp.run(f'cp {Mfit_dir}/exclude_w.dat {self.mfit_object_path}/exclude_w.dat', shell=True)
        logger.debug(f'Command {Command3.args} ran {"successfully" if Command3.returncode ==0 else "unsuccessfully"}')
        # logger.debug(f'Command {Command3.args} ran with returncode {Command3.returncode}')
        Command4 = sp.run(f'cp {Mfit_dir}/exclude_p.dat {self.mfit_object_path}/exclude_p.dat', shell=True)
        logger.debug(f'Command {Command4.args} ran {"successfully" if Command4.returncode ==0 else "unsuccessfully"}')
        # logger.debug(f'Command {Command4.args} ran with returncode {Command4.returncode}')
        logger.info('Molecfit parameter files prepared')


    def modify_par_file(self):
    #function to modify the template parameter file to add correct paths and files
        with self.par_path.open(mode='w') as fout:
            with Mfit_dir.joinpath('template.par').open(mode="r") as f:
                lines=f.readlines()
                for line in lines:
                    if line.strip().startswith('user_workdir'):
                        self.mfit_object_path
                        line = f'user_workdir: {self.mfit_object_path}'
                        logger.debug(f'user_workdir updated: {self.mfit_object_path}')
                    elif line.strip().startswith('filename'):
                        line=f'filename: {self.mfit_object_path}/00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits'
                        logger.debug(f'filename updated: {self.mfit_object_path}/00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits')
                    elif line.strip().startswith('wrange_include'):
                        line=f'wrange_include: {self.mfit_object_path}/include.dat'
                        logger.debug(f'wrange_include updated: {self.mfit_object_path}/include.dat')
                    elif line.strip().startswith('wrange_exclude'):
                        line=f'wrange_exclude: {self.mfit_object_path}/exclude_w.dat'
                        logger.debug(f'wrange_exclude updated: {self.mfit_object_path}/exclude_w.dat')
                    elif line.strip().startswith('prange_exclude'):
                        line=f'prange_exclude: {self.mfit_object_path}/exclude_p.dat'
                        logger.debug(f'pixel ranges updated: {self.mfit_object_path}/exclude_p.dat')
                    elif line.strip().startswith('output_dir'):
                        line=f'output_dir: {self.mfit_object_path}/output'
                        logger.debug(f'output_dir updated: {self.mfit_object_path}/output')
                    elif line.strip().startswith('output_name'):
                        line=f'output_name: {self.object}_{self.unseq}'
                        logger.debug(f'output name updated: {self.object}_{self.unseq}')
                    fout.write(line)

    def run_molecfit(self):
        clean = sp.run('rm -rf /tmp/mol*',shell=True)
        logger.info(f'{clean} ran with returncode {clean.returncode}')
        activate_comm = sp.run(f'source /home/karansinghd/.bashrc; conda activate py27; ~/Applications/Molecfit/bin/molecfit {self.par_path};/home/karansinghd/Applications/Molecfit/bin/calctrans {self.par_path};conda deactivate',shell=True)
        logger.info(activate_comm)

    def rewrite_fits(self):
        wave, flux, var, hd, wl_org = self.load_HERMES_data()
        # header = fits.getheader(f'/STER/mercator/hermes/{self.night}/reduced/00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_cf.fits')
        logger.debug('Clipping negative flux values')
        # #Correct for atm extinction
        real_flux = (flux>=0)
        logger.debug(f'{real_flux.sum()} number of pixels with positive fluxes are remaining')
        wave = wave[real_flux]
        flux = flux[real_flux]
        var = var[real_flux]
        wl_org = wl_org[real_flux]
        logger.debug('Correcting for atm extinction')
        ext_wl, aext_coeff = np.loadtxt('/STER/karansinghd/PhD/ResponseCorrection/V1/atm_ext.dat', dtype=float, unpack=True)
        ref_atm_cor = self.get_atm_ext_cor(ext_wl, aext_coeff, hd, wave)
        atm_cor_flux = flux * ref_atm_cor
        header=self.get_mfit_header(hd)
        #Write out new fits
        self.write_mfit_fits(wave, atm_cor_flux, var, header, wl_org)
        # exit()

    def get_mfit_header(self,hd):
        try:
            hd['MJD-OBS']
            logger.debug('Header already Molecfit compliant')
            return
        except:
            logger.debug('Header has to be edited!')
            hd['EXTNAME'] = ('FLUX','')
            hd['CD1_1']   = (float(hd["CDELT1"]),'Wavelength step taken from CDELT1')
            jd            = Time(hd['DATE-OBS'], format = 'isot', scale='utc').jd
            hd['MJD-OBS'] = (jd-2400000.5,'Modified Julian date (JD - 2400000.5)')
            utc           = np.array(hd['DATE-OBS'].split("T")[1].split(":")).astype(float)
            utc           = utc[0]*3600.+utc[1]*60.+utc[2]
            hd['UTC']     = (utc,'UTC start time of the exposure, in seconds')
            hd['TELALT']  = float(hd["TELALT"])
            hd['TEMPM1']  = float(hd["TEMPM1"])
            hd['HISTORY'] = "20180216: UTC,MJD-OBS,CD1_1,EXTNAME added for use in MOLECFIT "
            hd['HISTORY'] = "20180216: TELALT,TEMPM1 conversion to float for use in MOLECFIT "
            return hd

    def write_fits(self):
        '''write the final output fits'''
        # logger.info(f'Trying to open {self.filename}')
        hdus = fits.open(self.mfit_object_path / Path(f'00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits'))
        # head = self.header
        # head = hdus[0].header
        # flx = fits.getdata(self.mfit_object_path / Path(f'00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_cf.fits'))
        # head = fits.getheader(self.mfit_object_path / Path(f'00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits'))
        # hdus[0].header = head
        # hdus[0].data = flx
        hdus[0].header = self.header
        hdus[0].header['STDNIGHT'] = self.STD_NIGHT
        logger.info(f'Header updated for STDNIGHT = {self.STD_NIGHT}')
        hdus[0].header['HISTORY'] = f'spectrum corrected for the instrumental response function fit with a spline (K. Dsilva)'
        hdus[0].header['HISTORY'] = f'Keyword STDNIGHT included to indicate standard star observed on the same night (K.Dsilva)'

        #remove the primary image
        hdus[0].data=[]
        #remove the BINTABLE from Molecfit
        logger.debug(hdus.info())
        _bt = hdus.pop(1)
        #define and add new BINTABLE
        col_wave = fits.Column(name='wave', format='E',array=self.wave)
        col_flux = fits.Column(name='flux', format='E',array=self.rc_flux_no_TAC)
        col_err = fits.Column(name='err', format='E',array=self.rc_err_no_TAC)
        col_flux_molec = fits.Column(name='flux_molec', format='E',array=self.rc_flux)
        col_err_molec = fits.Column(name='err_molec', format='E',array=self.rc_err)

        cols = fits.ColDefs([col_wave,col_flux,col_err,col_flux_molec,col_err_molec])
        hdu_new = fits.BinTableHDU.from_columns(cols,name ='CRF_TABLE')
        hdus.append(hdu_new)
        hdus[1].header['HISTORY'] = 'Flatfield and instrumental response corrected spectrum saved in columns: wave, flux, err, flux_molec, err_molec'
        # logger.info(hdus.info())
        logger.info(f'New generated fits: {self.outpath}')
        hdus.writeto(self.outpath, overwrite=True)
        logger.info('Done!')

    def write_mfit_fits(self, wl, flx, var, head, wl_original):
        '''Writes a fits file that is molecfit compliant
        '''
        # head['EXTEND'] = True
        err = np.sqrt(var)
        ifn = Path(f'/STER/mercator/hermes/{self.night}/reduced/00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits')
        try:
            hdu = fits.open(ifn)
        except:
            logger.critical("File does not exist in Pierre's directory")
        hdu[0].header = head
        hdu[0].header['HISTORY'] = 'Flux corrected for atmospheric extinction (K. Dsilva)'
        hdu[0].data = flx
        #Create a new BINTABLE and remove the old
        col_wave = fits.Column(name='Wavelength', format='E',array=wl)
        col_flux = fits.Column(name='Flux', format='E',array=flx)
        col_err = fits.Column(name='Flux_err', format='E',array=err)
        col_wave_original = fits.Column(name='Wavelength_original', format='E', array=wl_original)
        cols = fits.ColDefs([col_wave_original,col_wave,col_flux,col_err])
        hdu_new = fits.BinTableHDU.from_columns(cols,name ='BINTABLE')
        while len(hdu)>1:
            a=hdu.pop()
        hdu.append(hdu_new)
        # __b = hdu.pop(1)
        # print(hdu.info())
        # exit()
        # hdu.insert(1,hdu_new)
        # logger.debug(hdu.info())
        hdu[1].header['HISTORY'] = 'BINTABLE added for molecfit usage'
        hdu[1].header['HISTORY'] = 'COLUMNS: Wavelength (w/o Bar Cor), Flux, Error and original Wavelength array (K. Dsilva)'
        # pl.plot(wl,flx)
        # pl.plot(wl_original,flx,'k')
        # pl.show()
        ofn = self.mfit_object_path / Path(ifn.stem + ifn.suffix)
        logger.debug(f'New generated fits: {ofn}')
        hdu.writeto(ofn, overwrite=True)

    def create_specific_plot(self,xdata,ydata,**kwargs):
        color=kwargs.pop('color','blue')
        linestyle=kwargs.pop('linestyle','-')
        linewidth=kwargs.pop('linewidth',3)
        xlabel=kwargs.pop('xlabel','X data')
        ylabel=kwargs.pop('ylabel','Y data')
        title=kwargs.pop('title','Figure1')
        #use the given stuff to set rcparams
        mpl.rc('lines',lw=linewidth,ls=linestyle,c=color)
        fig,ax1=pl.subplots()
        ax1.plot(xdata,ydata)
        # minimum = 0.9 * min(ydata[(xdata > 5000) & (xdata < 6000)])
        ax1.xaxis.set_tick_params(labelsize=12)
        ax1.yaxis.set_tick_params(labelsize=12)
        # ax1.set_ylim([0.8,10])
        ax1.set_xlabel(xlabel,fontsize=14)
        ax1.set_ylabel(ylabel,fontsize=14)
        # ax1.set_ylim((minimum,1))
        ax1.set_title(title)
        # pl.savefig(f'/STER/karansinghd/PhD/ResponseCorrection/testing/Plots/{title}.pdf',format='pdf')
        pl.savefig(f'{LocalPath}/{self.object}/{title}.pdf',format='pdf')
        pl.close()

    def create_plots(self,Resp,Poly):
        '''Function to create a plot of the response'''
        logger.info('Creating a plot of the response')
        mpl.rc('lines',lw=3)
        fig,(ax1,ax2)=pl.subplots(nrows=2, ncols=1, sharex=True, sharey=False, gridspec_kw={'height_ratios':[3,1]},figsize=(8,8))
        ax1.plot(self.wave,Resp,'r',label='Median filtered response')
        ax1.plot(self.wave,Poly,'k--',label=f'Spline fit')
        ax2.plot(self.wave,Resp-Poly,'k')
        ax2.axhline(y=0,color='r')
        ax2.xaxis.set_tick_params(labelsize=12)
        ax1.yaxis.set_tick_params(labelsize=12)
        ax2.yaxis.set_tick_params(labelsize=12)
        ax2.set_xlabel('Wavelength ($\AA$)',fontsize=14)
        ax1.set_ylabel('ADU/(erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$)',fontsize=14)
        ax2.set_ylabel('Residuals',fontsize=14)
        ax1.set_title(f'Response curve for NIGHT: {self.night}, OBJECT:{self.object}, UNSEQ:{self.unseq}')
        ax1.legend(loc='upper left')
        fig.tight_layout()
        Path('plots').mkdir(exist_ok=True)
        pl.savefig(f'plots/response_{self.night}_{self.unseq}.png')
        pl.close()

    def fit_spline(self,r1):
        x=self.wave
        y=r1

        left = max(4000,min(self.wave))
        right = min(9000,max(self.wave))
        #create knot array
        violet = np.linspace(left,4580,10)
        blue = np.linspace(4600,5525,25)
        green = np.linspace(5530,6635,30)
        red = np.linspace(6650,right,20)

        knots = np.concatenate((violet,blue,green,red))

        spl = us(x,y,t=knots)
        return(spl(x))

if __name__=='__main__':
    # correct_spectrum(20100930,307479,'V__AG_Psc')
    # correct_spectrum(20100927,307233,'NSV___118')
    # f=correct_spectrum(20111009,376039,'HIP_113222')
    # f = correct_spectrum(20110223,334866,'HD14055')
    f = correct_spectrum(20101104,314318,'NSV___691')
    exit()
    # df = pd.read_csv('prog28_good_20190702.csv',header=0,sep=',')
    df = pd.read_csv('melchiors_meta.csv',sep='|')
    # df = pd.read_csv('OverviewProgram28.csv',header=0,sep='|')
    # df.columns = [x.strip() for x in df.columns]
    # print('NEW',df.columns)
    '''
    snippet for specific numbers
    # '''
    # # nums_to_check =[]
    # for p in Path('/STER/pierre/hermes/p28/funkyblue/').glob('*.png'):
    #     unseq = int(p.stem.split('_')[0])
    #     nums_to_check.append(unseq)
        # print(unseq)
        # exit()
    # nums_to_check = [453599,453600,453601,453602,453603]
    # nums_to_check = [453585,453591,453593,453595,453596,453597]

    # ser = pd.read_csv('/STER/karansinghd/List_to_check.csv')
    # print(np.array(ser))
    # nums_to_check =np.array(ser)
    # exit()
    with Path('./log.log').open('w') as logfile:
        for i,row in df.iterrows():
            night = row['night']
            # night = int(row['filename'].split("/")[4])
            # object = " ".join(row['  starname  '].strip()).replace('*','_')
            object = row['starname'].strip().replace(' ','_').replace('*','_')
            unseq = int(row['obsid'])
            # if unseq not in nums_to_check:
            # if unseq not in [412343]:
                # continue
            # if (night < 20120101):
            #     logfile.write(f'{unseq},{night} - no meteo\n')
            #     # print(,)
                # continue
            # if night > 20150100:
                # exit()
            '''
            snippet to copy spectra
            '''
            # print(object,unseq,night)
            # shutil.copy(f'/STER/karansinghd/PhD/Projects/P28_c/{object}/corrected/00{unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits',f'/STER/karansinghd/P28_c_cr/00{unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits')
            # print('File copied')
            # continue
            try:
                print('-'*73)
                f=correct_spectrum(night,unseq,object)
                df_stdinfo.loc[i]=[unseq,f.stdunseq,f.stdname,night,f.stdnight]
                shutil.copy(f'/STER/karansinghd/PhD/Projects/P28_c/{object}/corrected/00{unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits',f'/STER/karansinghd/P28_c_cr/00{unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_cr.fits')
                print('-'*73)
                # exit()
            except Exception as e:
                logfile.write(f'{logger.critical(traceback.format_exc())}\n Error: {e}\n')
                print(f'Error: {e}')
                # exit()
    df_stdinfo.to_csv('stdinfo.csv',index=False)
