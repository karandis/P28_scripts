'''
Hermes Response Correction v1.0
Author : Karan Dsilva (karan.dsilva@kuleuven.be)
UPDATE: 28/01/2022:
        mean RV used to correct for all P2 stars except GSC-

METHOD : Takes in a night as input and check:

    a. if there is a Program 2 star observed that night
    b. if yes, is it a star for which we have a good model
    c. if no to any of those questions, look for the closest P2 star in time (within the provided tolerance)
    d. determine the instrumental response for the night
    e. save it to text and plot the response

'''
#------------------------------- IMPORTS  --------------------------------------
import datetime as dt
from pathlib import Path
import numpy as np
import scipy.interpolate as si
import pandas as pd
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as pl
import plotting
plotting.set_mode('paper')
import scipy.constants as sc
from astropy.io import fits
from scipy.signal import medfilt as mf
from numpy.polynomial import polynomial as poly
import subprocess as sp
import os
import logging
from astropy.time import Time
import RegscorePy as rp #Module to get the Bayesian Information Criterion
from scipy.interpolate import LSQUnivariateSpline as us
import traceback
#---------------------------- REFERENCE SEDs -----------------------------------
#The reference SEDs are corrected for reddening (interstellar extinction)
#scripts to redden the SEDs can be found in /STER/karansinghd/PhD/ResponseCorrection/V1/InputFiles/
ref_sed = {'HD152614': 'HD152614.rgs',                                       # B
           'HD36267': 'HD36267.rgs', 'hd36267': 'HD36267.rgs',       # B
           'HD118098': 'HD118098.rgs', 'HD 118098': 'HD118098.rgs',  # A
           'HD14055': 'HD14055.rgs', 'hd14055': 'HD14055.rgs',       # A
           'HD87887': 'HD87887.rgs',                                         # A
           'HD46300': 'HD46300.rgs', 'hd46300': 'HD46300.rgs', 'HD 46300': 'HD46300.rgs', # A
           'HD184006': 'HD184006.rgs',                                       # A
           'HD149212': 'HD149212.rgs',                                       # A
           'HD147449': 'HD147449.rgs',                                       # A
           'GSC4293-0432': 'GSC4293-0432.rgs',                               # A
           'HD214994': 'HD214994.rgs',                                       # A
           'HD42818': 'HD42818.rgs', 'HD 42818': 'HD42818.rgs',  # A
           'HD56169': 'HD56169.rgs',                                         # A
           # 'HD185395': 'HD185395.rgs',                                       # F
           # 'HD206826': 'HD206826.rgs',                                       # F
           # 'HD220657': 'HD220657.rgs',                                       # F
           'HD84937': 'HD84937.rgs',
           #'BD+17.4708': 'BD+174708.rgs'                                     # F
           }


RVs  = {'HD14055':6.586,
        'HD214994':9.420,
        'HD185395':-27.369,
        'HD152614':-16.725,
        'HD220657':-10.029,
        'HD46300':12.752,
        'HD184006':-14.398,
        'GSC4293-0432':-3.491,
        'HD56169':-2.591,
        'HD118098':-10.086,
        'HD84937':-14.920,
        'HD149212':-6.338,
        'HD36267':	21.487,
        'HD206826':18.227,
        'HD87887':7.023,
        'HD42818':-3.899,
        'HD147449':-45.706}
#----------------------------------------------------------------
logging.basicConfig(level=logging.INFO,format='%(asctime)s: [%(name)s:%(levelname)s] %(message)s',datefmt='%d-%m %H:%M')
logger = logging.getLogger("RESPONSE")
#----------------------------------------------------------------
class response:
    '''Uses Molecfit corrected fits files
    this class has the following variables, which can be accessed with self.* :
        #Initialization
        Night: main night being processed
        tolerance: number of days to look in time
        night: the night of the response being calculated (can change based on
               tolerance)
        unseq: for the unique observation

        #I/O
        HermesPath: Path to the HERMES database
        Mfit_dir: Path to the main Molecfit directory
        InputFilesDir: Path to model SEDs
        rpath:  Path to the saved response
        filename: Path to the Molecfit corrected spectrum

        #Molecfit related
        wavelength_original_index: index for the wl array output from Molecfit
        flux_TAC_index: index for the Molecfit corrected flux
        flux_index: index for the original flux (used for RV calculation)
        mfit_object_path: Path to the Molecfit folder containing
                          observations of that specific star
        par_path: Path to the Molecfit input file
        redo_molec: keyword to redo the Molecfit run

        #Fallbacks and indicators
        flag: to indicate if the spectrum is reduced or if it is only in the HERMES Overview file (problem with reduction).
        key: 1 if meteo data is available else 0
        list_failed: list of failed reductions for that night given tolerance

        #The data
        wave: the original wavelength array (NOTE: negative  flux pixels clipped)
        flux: the original (uncorrected) flux array
        flux_corr: the telluric absorption corrected (TAC) flux
        rv: Radial velocity (RV) from dictionary calculated by Hans (2022 Jan)
        model_flux: flux from the model SED
        response: the median-filtered response (flux_corr/model_flux)

        #Spline fitting stuff
        knots: knotpoints used to fit the spline for the response.
        spline_fit: the spline fit to the response (NEW updated)
        residuals: residuals of the fit

    '''
    #Indices for load_FITS_table_data
    wavelength_original_index = 0
    flux_TAC_index = 5
    flux_index = 2
    #Paths to frequently used directories
    HermesPath = Path('/STER/mercator/hermes')
    Mfit_dir = Path('/STER/karansinghd/PhD/ResponseCorrection/Molecfit')
    InputFilesDir = Path('/STER/karansinghd/PhD/ResponseCorrection/ModelSEDs')

    def __init__(self,night,tolerance=None,overwrite=False):
        self.redo_molec=False
        self.tolerance = tolerance
        self.Night = night
        self.list_failed = []
        print('--------------------------------------------------')
        logger.info(f'Processing NIGHT: {self.Night}')
        print('--------------------------------------------------')
        p2_info= self.get_p2_info()
        p2_info.reset_index(inplace=True,drop=True)
        if not len(p2_info):
            logger.critical('No Program 2 star for the night!')
            return
        for j,row in p2_info.iterrows():
            global i
            i = j
            row = p2_info.iloc[i]
            self.unseq = row['unseq']
            #remove one bad model
            if self.unseq == 325226:
                continue
            self.object = row['object']
            self.night = row['night']
            self.rpath = Path(f'/STER/karansinghd/PhD/ResponseCorrection/responses_c_2022/{self.night}_{self.object}_{self.unseq}.txt')
            print('----------------------------------------------------')
            logger.info(f'NIGHT: {self.night}, UNSEQ: {self.unseq}, OBJECT: {self.object}')
            self.mfit_object_path = self.Mfit_dir.joinpath(self.object)
            self.par_path = self.mfit_object_path/Path(f'{self.object}_{self.unseq}.par')
            self.filename = Path(f'{self.mfit_object_path}/output/00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c_TAC.fits')
            #Run Molecfit if not already done
            if (not self.filename.is_file()) or (self.filename.is_file() and self.redo_molec):
                logger.info('Molecfit will be run as TAC file does not exist')
                self.prepare_molecfit_files()
                self.modify_par_file()
                try:
                    self.rewrite_fits()
                    if not self.key:
                        logger.critical('Keywords for meteo station readings missing in header, next!')
                        continue
                except Exception as e:
                    logger.warning('Failed while rewriting FITS file')
                    if not (self.flag and self.flag_var):
                        e = f'One of the input files is missing - flux: {self.flag} and/or variance: {self.flag_var}'
                    logger.critical(f'Exception: {e}')
                    self.list_failed.append([self.night,self.unseq,e])
                    continue
                self.run_molecfit()
            try:
                self.wave, self.flux, self.flux_corr = self.load_FITS_table_data()
                logger.info('Successfully loaded Molecfit corrected spectrum')
            except Exception as e:
                logger.critical(f'Problem with the Molecfit corrected spectrum: {e}')
                self.list_failed.append([self.night,self.unseq,e])
                continue
            self.rv = self.calc_rad_vel()
            self.wave = self.correct_vel_shift()
            self.model_flux = self.load_model_SED()
            if self.rpath.is_file() and not overwrite:
                logger.info(f'{self.rpath} exists')
                data = pd.read_csv(self.rpath,sep='\t',header=0,encoding='utf-8',engine='python')
                self.wave = data.wavelength.to_numpy()
                self.response = data.response.to_numpy()
                self.spline_fit = data.spline.to_numpy()
                self.residuals =  self.response - self.spline_fit
            else:
                self.response = self.get_response()
                self.spline_fit = self.fit_spline()
                self.residuals =  self.response - self.spline_fit
                logger.info(f'Spline fit was successful.')
                self.create_plots()
                self.save_response()
                logger.info(f'Response successfully saved to {self.rpath}')
            return


    def load_FITS_table_data(self):
        '''
        Function to load table data from the FITS file output by Molecfit.
        '''
        hdul = fits.open(self.filename)
        head = fits.getheader(self.filename)
        table_data = hdul[1].data
        wave = table_data.field(self.wavelength_original_index)
        flux = table_data.field(self.flux_index)
        flux_TAC = table_data.field(self.flux_TAC_index)
        return wave, flux, flux_TAC

    def calc_rad_vel(self):
        '''Function to calculate the RV of the spectrum using a CCF
        Now that we drop BD+17, we can reduce the  CCF exploration velocity from
        +/-320 Km/s to 150 Km/s
        The Balmer lines are used to derive the RV (not ideal but we do not have
        a normalized spectrum to do a better CCF)
        We do not use the Gaussian fitting but the bisector analysis, which is
        much more reliable

        Update 28/01/2022: Scrapping the majority of this function. No need to calc RVs, we just read them from a file!
        '''
        global i
        if self.object == 'GSC4293-0432':
            #read RV from file!
            rvdat =  pd.read_csv('/STER/karansinghd/PhD/ResponseCorrection/RadialVelocitiesGSC4293-0432.vrdata',delim_whitespace=True,header=None)
            rvdat.columns=['JD','vrad','sigma','unseq']
            rv = rvdat.loc[rvdat.unseq==self.unseq].vrad.to_numpy()[0]
        else:
            rv = RVs[self.object]

        return rv

    def correct_vel_shift(self):
        '''Function to correct for the RV shift before getting  the  response'''
        radvel = self.rv*1000 #in m/s
        logger.info(f'Correcting for radial velocity of:{self.rv:.2f} Km/s')
        return self.wave*(1+radvel/sc.c)

    def load_model_SED(self):
        '''Function to load the model SED of the corresponding program 2 star'''
        logger.debug(f'Model Path: {self.InputFilesDir}/{ref_sed[self.object]}')
        ref_wave, ref_flux = np.genfromtxt(fname=f'{self.InputFilesDir}/{ref_sed[self.object]}', unpack=True)
        #interpolate the model at the wavelengths of the spectrum
        model_flux = si.interp1d(ref_wave, ref_flux, kind='linear',fill_value='extrapolate')(self.wave)
        logger.info('Successfully loaded model SED')
        return model_flux

    def get_response(self):
        '''Function to get the response and apply a median filter '''
        self.rough_response = self.flux_corr/self.model_flux
        #2 median filters to deal with the edges
        R1_1 = mf(self.rough_response,1001)
        R1_2 = mf(self.rough_response,11)
        R1_1[:200] = R1_2[:200]
        R1_1[-40:] = R1_2[-40:]
        return R1_1


    def fit_spline(self):
        x=self.wave
        y=self.response

        left = min(x)+10
        right = max(x)-10

        #create knot array
        added = np.array([3781,3852,3913,3952])
        if self.object =='HD84937':
            violet = np.linspace(4000,4340,25)
        else:
            violet = np.linspace(4000,4280,20)

        blue = np.linspace(4350,5515,30)
        green = np.linspace(5545,6635,30)
        red = np.linspace(6675,right,20)

        knots = np.concatenate((added,violet,blue,green,red))
        self.knots=knots

        spl = us(x,y,t=knots)
        return(spl(x))

    def save_response(self):
        '''Function to save the response as a text file'''
        savefile=pd.DataFrame({'wavelength':self.wave,'spline':self.spline_fit,'response':self.response})
        savefile.to_csv(self.rpath,header=True,index=False,sep='\t')

    def prepare_molecfit_files(self):
        '''Function to create Molecfit input directories and input files.
        Independent directories are created in order to offer traceback and
        modification of the parameter file to run Molecfit again
        '''
        if not os.path.isdir(self.mfit_object_path):
            os.mkdir(self.mfit_object_path)
        if not os.path.isdir(self.mfit_object_path/Path('output')):
            os.mkdir(self.mfit_object_path/Path('output'))
        #copy necessary parameter files
        command1 = sp.run(f'cp {self.Mfit_dir}/include.dat {self.mfit_object_path}/include.dat', shell=True)
        logger.debug(f'Command {command1.args} ran with returncode {command1.returncode}')
        command2 = sp.run(f'cp {self.Mfit_dir}/exclude_w.dat {self.mfit_object_path}/exclude_w.dat', shell=True)
        logger.debug(f'Command {command2.args} ran with returncode {command2.returncode}')
        command3 = sp.run(f'cp {self.Mfit_dir}/exclude_p.dat {self.mfit_object_path}/exclude_p.dat', shell=True)
        logger.debug(f'Command {command3.args} ran with returncode {command3.returncode}')
        logger.info('Molecfit inclusion and exclusion files prepared')


    def modify_par_file(self):
        '''Function to modify the template parameter file to add correct paths and files
        '''
        with self.par_path.open(mode='w') as fout:
            with self.Mfit_dir.joinpath('template.par').open(mode="r") as f:
                lines=f.readlines()
                for line in lines:
                    if line.strip().startswith('user_workdir'):
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
        logger.info("Molecfit parameter file created successfully")

    def run_molecfit(self):
        '''Function to run molecfit. As it runs on  Python 2, a new environment
        has to be activated. It creates its workspace in the tmp directory, and
        hence has to be cleaned in order to avoid the code crashing
        '''
        #Command to remove any previous molecfit runs
        clean = sp.run('rm -rf /tmp/mol*',shell=True)
        logger.debug(f'{clean} ran with returncode {clean.returncode}')
        activate_comm = sp.run(f'source /home/karansinghd/.bashrc; conda activate py27; /home/karansinghd/Applications/Molecfit/bin/molecfit {self.par_path};/home/karansinghd/Applications/Molecfit/bin/calctrans {self.par_path};conda deactivate',shell=True)
        logger.info(f'Molecfit ran with a returncode {activate_comm.returncode} (0 = successful)')
        logger.debug(activate_comm)

    def load_HERMES_data(self):
        ''' Function to load a HERMES spectrum.
        Returns:
            wavelength : numpy array of wavelengths(lambda or ln(lambda)) in Angstrom
            flux : numpy array with flux (photons?)
            var : variance of the flux
            head : header of the fits file
            Wavelength_original: barycentric correction included
        '''
        #The input file names
        ifn = self.HermesPath.joinpath(f'{self.night}/reduced/00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits')
        ifn_old = self.HermesPath.joinpath(f'{self.night}/reduced/{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits')

        ifn_var = self.HermesPath.joinpath(f'{self.night}/reduced/00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_mergedVar_c.fits')
        ifn_var_old = self.HermesPath.joinpath(f'{self.night}/reduced/{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_mergedVar_c.fits')
            #The data
        try:
            flux, head = fits.getdata(ifn,header=True)
            self.flag=1
        except:
            try:
                flux, head = fits.getdata(ifn_old,header=True)
                self.flag=1
            except:
                self.flag=0
                logger.critical(f'Cannot find the spectrum (UNSEQ:{self.unseq}, NIGHT:{self.night}) in the HERMES database!')
        try:
            var = fits.getdata(ifn_var)
            self.flag_var=1
        except:
            try:
                var = fits.getdata(ifn_var_old)
                self.flag_var=1
            except:
                self.flag_var=0
                logger.critical(f'Cannot find the variance (UNSEQ:{self.unseq}, NIGHT:{self.night}) in the HERMES database!')
        try:
            head['HUM_MET']
            self.key=1
        except KeyError:
            self.key=0
            logger.warning('FITS header does not contain Meteo station keywords')
        logger.debug('FITS files loaded, condensing data to arrays')
        crpix = head['CRPIX1']-1
        n_points = head['NAXIS1']
        delta_w = head['CDELT1']
        start_w = head['CRVAL1'] - crpix * delta_w
        offset = self.correct_bar_vel_shift(head)
        start_w_corr = start_w+ offset
        logger.info('Barycentric velocity correction has been undone')
        #create both the original and the corrected wavelengths to be saved in the BINTABLE
        log_wl_original = np.linspace(start_w, start_w+(delta_w * (n_points-1)), n_points)
        log_wl = np.linspace(start_w_corr, start_w_corr+(delta_w * (n_points-1)), n_points)
        #We need them in linear scale
        wave = np.exp(log_wl)
        wave_original = np.exp(log_wl_original)
        #Check for and clean 1) nans and 2) zeros (ORDER IS IMPORTANT:  T.Merle)
        logger.debug('Cleaning nans and zeros')
        wave, flux, var, wave_original = self.clean_nan(wave,flux,var,wave_original)
        logger.debug('nans cleaned')
        wave, flux, var, wave_original = self.clean_zeros(wave,flux,var,wave_original)
        logger.debug('zeros cleaned')
        return wave, flux, var, head, wave_original

    def write_fits(self, wl, flx, var, head, wl_original):
        '''Called by rewrite_fits() to save a FITS file that is molecfit compliant
        '''
        err = np.sqrt(var)
        ifn = self.HermesPath.joinpath(f'{self.night}/reduced/00{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits')
        ifn_old = self.HermesPath.joinpath(f'{self.night}/reduced/{self.unseq}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits')
        try:
            hdul = fits.open(ifn)
        except:
            try:
                hdul = fits.open(ifn_old)
            except:
                logger.critical(f'Cannot find the spectrum (UNSEQ:{self.unseq}, NIGHT:{self.night}) in the HERMES database!')
                return
        hdul[0].header = head
        hdul[0].data = flx
        #Create a new BINTABLE
        col_wave = fits.Column(name='Wavelength', format='E',array=wl)
        col_flux = fits.Column(name='Flux', format='E',array=flx)
        col_err = fits.Column(name='Flux_err', format='E',array=err)
        col_wave_original = fits.Column(name='Wavelength_original', format='E', array=wl_original)
        cols = fits.ColDefs([col_wave_original,col_wave,col_flux,col_err])
        hdu_new = fits.BinTableHDU.from_columns(cols,name ='BINTABLE')
        logger.debug('Created the new HDU from the columns')
        #in case the InstrumentConfig extension is still there, we need to remove it!
        while len(hdul)>1:
            a=hdul.pop()
            logger.debug("Removed the InstrumentConfig extension")
        hdul.append(hdu_new)
        logger.debug('Appended the new HDU')
        hdul[1].header['HISTORY'] = 'Flux corrected for atmospheric extinction (K. Dsilva)'
        hdul[1].header['HISTORY'] = 'BINTABLE added for molecfit usage'
        hdul[1].header['HISTORY'] = 'COLUMNS: Wavelength (w/o Bar Cor), Flux, Error and original Wavelength array (K. Dsilva)'
        ofn = self.mfit_object_path / Path(ifn.stem + ifn.suffix)
        logger.debug(f'New generated fits: {ofn}')
        hdul.writeto(ofn, overwrite=True)
        logger.info('Molecfit compliant FITS file written')

    def get_atm_ext_cor(self,ext_wl, aext_coeff, head, wl):
        '''Function to get the atmospheric extinction (adapted from T.Merle)
        NOTE: It is airmass dependent
        '''
        m = (min(wl) < ext_wl) & (ext_wl < max(wl))
        ref_aext = si.interp1d(ext_wl[m], aext_coeff[m], kind='quadratic', bounds_error=False, fill_value='extrapolate')
        z = np.deg2rad(90 - float(head['TELALT'])) # Zenithal angle in radians
        logger.debug(f'Zenithal distance: {z:2.5f} rad')
        atm_cor = 10**(0.4*ref_aext(wl)/np.cos(z))
        return atm_cor

    def clean_nan(self,x0, y0, z0, x1):
        '''Function to clean NaNs (adapted from T.Merle)
        '''
        m = (y0 != y0) | (z0 != z0)

        x = x0[~m]
        y = y0[~m]
        z = z0[~m]
        x2 = x1[~m]

        try:
            y = si.interp1d(x, y, kind='linear', bounds_error=True)(x0)
            z = si.interp1d(x, z, kind='linear', bounds_error=True)(x0)
            x = x0
            x2 = x1
        except ValueError:
            pass

        return x, y, z, x2

    def clean_zeros(self,x0, y0, z0, x1):
        '''Function to clean zeros (adapted from T.Merle)
        '''
        m = (y0 == 0) | (z0 == 0)

        x = x0[~m]
        y = y0[~m]
        z = z0[~m]
        x2 = x1[~m]

        y = si.interp1d(x, y, kind='linear', fill_value='extrapolate')(x0)
        z = si.interp1d(x, z, kind='linear', fill_value='extrapolate')(x0)

        return x0, y, z, x1

    def correct_bar_vel_shift(self,Head):
        '''Function to undo barycentric correction.
        Input:
            Header of fits (to get BVCOR): Km/s
        Returns:
            offset to correct the wavelength array
        '''
        if 'BVCOR' in Head:
            bar_v = Head['BVCOR'] #in Km/s
        else:
            logger.warning('BVCOR not found, assuming 0')
            bar_v = 0
        #Doppler shift corrected
        logger.info(f'Barycentric velocity from header: {bar_v} km/s')
        offset = np.log(1-(bar_v)/3.0e5)
        return offset


    def rewrite_fits(self):
        '''Function to correct for atmospheric extinction, remove pixels with
        negative flux, get a Molecfit compliant header and  write it out
        '''
        wave, flux, var, header, wl_org = self.load_HERMES_data()
        #Correct for atm extinction
        real_flux = (flux>=0)
        logger.debug(f'{real_flux.sum()} number of pixels are remaining')
        wave = wave[real_flux]
        flux = flux[real_flux]
        var = var[real_flux]
        wl_org = wl_org[real_flux]
        ext_wl, aext_coeff = np.loadtxt('/STER/karansinghd/PhD/ResponseCorrection/V1/atm_ext.dat', dtype=float, unpack=True)
        ref_atm_cor = self.get_atm_ext_cor(ext_wl, aext_coeff, header, wave)
        atm_cor_flux = flux * ref_atm_cor
        logger.debug('Corrected for atmospheric extinction')
        header=self.get_mfit_header(header)
        logger.debug('Header is compliant with Molecfit')
        #Write out new fits
        self.write_fits(wave, atm_cor_flux, var, header, wl_org)

    def get_mfit_header(self,hd):
        '''Function to generate a Molecfit compliant header (if it does  not exist)
        '''
        hd['EXTNAME'] = ('FLUX','')
        hd['CD1_1']   = (float(hd["CDELT1"]),'Wavelength step taken from CDELT1')
        jd            = Time(hd['DATE-OBS'], format = 'isot', scale='utc').jd
        hd['MJD-OBS'] = (jd-2400000.5,'Modified Julian date (JD - 2400000.5)')
        utc           = np.array(hd['DATE-OBS'].split("T")[1].split(":")).astype(float)
        utc           = utc[0]*3600.+utc[1]*60.+utc[2]
        hd['UTC']     = (utc,'UTC start time of the exposure, in seconds')
        hd['TELALT']  = float(hd["TELALT"])
        hd['TEMPM1']  = float(hd["TEMPM1"])
        hd['HISTORY'] = f'{dt.date.today().isoformat()}: UTC,MJD-OBS,CD1_1,EXTNAME added for use in MOLECFIT '
        hd['HISTORY'] = f'{dt.date.today().isoformat()}: TELALT,TEMPM1 conversion to float for use in MOLECFIT '
        return hd

    def get_p2_info(self):
        df = pd.read_csv("/STER/mercator/hermes/HermesFullDataOverview.tsv",sep='\t',skiprows=2)
        header = ['unseq', 'prog_id', 'obsmode', 'bvcor', 'observer', 'object', 'ra','dec', 'bjd', 'exptime', 'pmtotal', 'date-avg', 'airmass','filename']
        df.columns=header
        nights =[int(df['filename'].values[i].split('/')[4]) for i in range(len(df))]
        df=df.assign(night=nights)
        #get program 2 stars for the night
        df_night = df.loc[(df['prog_id']==2) & (df['night']==self.Night)]
        #Check if we have a good model of the p2 star(s)
        indices_good_models=self.check_obs_nights(df_night['object'])
        df_good_models=df_night.loc[indices_good_models]
        #If we don't, look for +/- tolerance
        if len(df_good_models)>0:
            logger.info(f'{len(df_good_models)} observations available for {self.Night}')
        else: #no good models, check if tolerance!
            if not self.tolerance:
                return df_good_models
            date = dt.date(year=int(self.Night/10000),month=int(self.Night/100)%100,day=self.Night%100)
            tol = dt.timedelta(days=self.tolerance)
            date_plus = (date+tol).isoformat()
            date_minus = (date-tol).isoformat()
            logger.info(f'No observation for which we have a good model. Looking for prog 2 stars between {date_plus} and {date_minus}')
            _s = ''
            date_plus = int(_s.join((date+tol).isoformat().split('-')))
            date_minus = int(_s.join((date-tol).isoformat().split('-')))
            df_tolerance = df.loc[(df['prog_id']==2) & ((df['night']<=date_plus) & (df['night']>=date_minus))]
            indices_good_models=self.check_obs_nights(df_tolerance['object'])
            df_good_models=df_tolerance.loc[indices_good_models]
        #If successful, we return the data frame with the info
            if len(df_good_models)>0:
                df_good_models.reset_index(drop=True,inplace=True)
                logger.info(f'{len(df_good_models)} observations available for {self.Night}, returning the closest')
                df_good_models['time_delta'] = (df_good_models['night']-self.Night).abs()
                df_good_models.sort_values(by='time_delta',ascending=True,inplace=True)
                df_good_models.filename = df_good_models.filename.astype(str)
            else:
                logger.critical(f'No good models available within +/- {self.tolerance} days of {self.Night}. Try again!')
        return df_good_models

    def check_obs_nights(self,series_objects):
        '''Function to check if the Prog2 stars are ones for which we have a
        good model of'''
        list_indices=[]
        for k,obj in series_objects.iteritems():
            if obj in ref_sed.keys():
                list_indices.append(k)
        return list_indices

    def create_plots(self):
        '''Function to create a plot of the response'''
        logger.info('Creating a plot of the fit')

        fig,(ax3,ax4,ax1,ax2)=pl.subplots(nrows=4, ncols=1, sharex=True, sharey=False, gridspec_kw={'height_ratios':[1,1,1,1]},figsize=(8,8))
        ax3.plot(self.wave,self.flux_corr)
        ax4.plot(self.wave,self.model_flux,zorder=12)
        ax1.plot(self.wave,self.rough_response,'b',label='Rough response',alpha=0.5)
        ax1.plot(self.wave,self.response,'r',label='Median filtered response')
        ax1.plot(self.wave,self.spline_fit,'k--',label=f'Spline fit')
        ax2.plot(self.wave,100*np.divide(self.residuals,self.spline_fit),'k')
        ax2.axhline(y=0,color='r')

        [ax1.axvline(x=j,color='grey',alpha=0.8) for j in self.knots]

        ax2.set_xlabel('Wavelength ($\AA$)')
        ax1.set_ylabel('ADU')
        ax2.set_ylabel('Residuals (%)')

        ax3.set_title(f'NIGHT: {self.night}, OBJECT:{self.object}, UNSEQ:{self.unseq}')
        fig.tight_layout()
        pl.savefig(f'/STER/karansinghd/PhD/Projects/P28_c/plots_2022/{self.night}_{self.object}_{self.unseq}.png',format='png')
        # pl.show()
        pl.close()
        logger.info(f'Figure saved to /STER/karansinghd/PhD/Projects/P28_c/plots_2022/{self.night}_{self.object}_{self.unseq}.png')


if __name__=='__main__':
    # 20100927_HD220657_307231
    #20110111_HD84937_327107
    #363584  HD152614  20110811.0  20110811.0
    #314382    307944  HD184006  20101105.0  20101005.0
    #374008    373993  HD185395  20110921.0  20110921.0
    #307630    307624  HD206826  20101002.0  20101002.0
    #389579    389569  GSC4293-0432  20111211.0  20111211.0
    # x=response(night=20111211,tolerance=None,overwrite=True)
    # x = response(night=20101220,tolerance=None,overwrite=True)
    # x = response(night=20111117,tolerance=60,overwrite=True)
    x = response(night=20120609,tolerance=None,overwrite=True)
    # correct_spectrum(20101220,325228,'HD36267')

    # HD36267, STDNIGHT: 20101220, STDUNUSEQ: 325226

    exit()
    '''
    df = pd.read_csv('melchiors_meta.csv',sep='|')
    df = df.drop_duplicates(subset=['night'])
    '''
    df = pd.read_csv('/STER/karansinghd/PhD/Projects/P28_c/stdinfo.csv')
    df.columns = ['unseq','stdunseq','stdname','night','stdnight']
    df = df.drop_duplicates(subset='night')
    df.reset_index(drop=True,inplace=True)
    for _,row in df.iterrows():
        print(f"--------------- ({_+1}/{len(df)}) ----------------")
        x = response(night=int(row['night']),tolerance=60,overwrite=False)
    print(df.head())
    exit()
    # with Path('./response.log').open('w') as logfile:
    #     startdate=dt.date(year=2010,month=1,day=1)
    #     delta = dt.timedelta(days=1)
    #     enddate = dt.date(year=2011,month=12,day=31)
    #     print(startdate,enddate)
    #     while startdate <= enddate:
    #         try:
    #             x=response(night=int(startdate.strftime(format="%Y%m%d")),tolerance=None,overwrite=True)
    #             startdate += delta
    #         except Exception as e:
    #             logfile.write(f'{startdate.strftime(format="%Y%m%d")},{logger.critical(traceback.format_exc())}\n Error: {e}\n')
    #             print(e)
    #             startdate += delta
