'''
Hermes Response Correction v1.0
Author : Karan Dsilva (karan.dsilva@kuleuven.be)
Adapted with input from T. Merle
UPDATE: *BD+174708 has been dropped due to poor corrections
        *The saved response has the polynomial as well as the median filtered
        saved.
        *saving the response with pandas, tab delimited
        *datetime corrected when looking for prog 2 stars using tolerance
        *normalization from  1 to 3 instead of -1 to 1
        *residuals plotted
        *Polynomial chosen based on Bayesian Information Criterion
        (BIC). Degree 30 has the lowest BIC and hence is the best model.
        *Both responses of the night kept in the form of arrays
        *added comments on all class variables

METHOD : Takes in a night as input and check:

    a. if there is a Program 2 star observed that night
    b. if yes, is it a star for which we have a good model
    c. if no to any of those questions, look for the closest P2 star in time
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
# mpl.use('agg')
import matplotlib.pyplot as pl
import plotting
plotting.set_mode('paper')
import scipy.constants as sc
from astropy.io import fits
from scipy.optimize import leastsq
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
           'HD185395': 'HD185395.rgs',                                       # F
           'HD206826': 'HD206826.rgs',                                       # F
           'HD220657': 'HD220657.rgs',                                       # F
           'HD84937': 'HD84937.rgs',
           #'BD+17.4708': 'BD+174708.rgs'                                     # F
           }

spec_types = {'HD152614': 'B',                                       # B
           'HD36267': 'B', 'hd36267': 'B',       # B
           'HD118098': 'A', 'HD 118098': 'A',  # A
           'HD14055': 'A', 'hd14055': 'A',       # A
           'HD87887': 'A',                                         # A
           'HD46300': 'A', 'hd46300': 'A', 'HD 46300': 'A', # A
           'HD184006': 'A',                                       # A
           'HD149212': 'A',                                       # A
           'HD147449': 'A',                                       # A
           'GSC4293-0432': 'A',                               # A
           'HD214994': 'A',                                       # A
           'HD42818': 'A', 'HD 42818': 'A',  # A
           'HD56169': 'A',                                         # A
           'HD185395': 'F',                                       # F
           'HD206826': 'F',                                       # F
           'HD220657': 'F',                                       # F
           'HD84937': 'F',
           #'BD+17.4708': 'BD+174708.rgs'                                     # F
           }
fname = {'A':'A_68Tau_depth0.95_blend20.list','B':'B_HR7512_depth0.98_blend20.list','F':"F_Procyon_depth0.95_blend20.list"}
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
        Mask_path: Path to list of Balmer lines for RV correction
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

        #Fallbacks and indicators
        key: 1 if meteo data is available else 0
        list_failed: list of failed attempts with unseq and error message

        #The data
        wave: the original wavelength array (NOTE: negative  flux pixels clipped)
        flux: the original (uncorrected) flux array
        flux_corr: the telluric absorption corrected (TAC) flux
        rv: Radial velocity (RV) from bisector analysis (CCF) using Balmer lines
        err_rv: the error on the RV
        model_flux: flux from the model SED
        response: the median-filtered response (flux_corr/model_flux)

        #Polynomial fitting stuff
        wl_norm: normalized wavelength array for polynomial fit
        spline_fit: the spline fit to the response (NEW updated)
        poly_degree: degree of the polynomial
        residuals: residuals of the fit
        rms: RMS value of the residuals
        bic: the Bayesian Information Criterion (see code for formula)


    '''
    #Indices for load_FITS_table_data
    wavelength_original_index = 0
    flux_TAC_index = 5
    flux_index = 2
    #Paths to frequently used directories
    HermesPath = Path('/STER/mercator/hermes')
    Mfit_dir = Path('/STER/karansinghd/PhD/ResponseCorrection/Molecfit')
    Mask_path = Path('/STER/karansinghd/PhD/ResponseCorrection/VelocityCorrection/masks/Balmer.list')
    InputFilesDir = Path('/STER/karansinghd/PhD/ResponseCorrection/ModelSEDs')
    LinesDir = Path('/STER/karansinghd/PhD/ResponseCorrection/VelocityCorrection/masks/')

    #Initialize the variables we want to use outside of the class as arrays
    night = dict()
    object = dict()
    unseq = dict()
    rpath = dict()
    par_path = dict()
    filename = dict()
    wave = dict()
    response = dict()
    residuals = dict()
    spline_fit = dict()
    wl_norm = dict()
    rv = dict()
    err_rv = dict()
    model_flux=dict()


    poly_degree = 30

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
            self.unseq[i] = row['unseq']
            self.object[i] = row['object']
            self.night[i] = row['night']
            self.rpath[i] = Path(f'/STER/karansinghd/PhD/ResponseCorrection/responses_c/{self.night[i]}_{self.object[i]}_{self.unseq[i]}.txt')
            if self.rpath[i].is_file() and not overwrite:
                logger.info(f'{self.rpath[i]} exists')
                continue
            print('##########################################################')
            logger.info(f'NIGHT: {self.night[i]}, UNSEQ: {self.unseq[i]}, OBJECT: {self.object[i]}')
            self.mfit_object_path = self.Mfit_dir.joinpath(self.object[i])
            self.par_path[i] = self.mfit_object_path/Path(f'{self.object[i]}_{self.unseq[i]}.par')
            self.filename[i] = Path(f'{self.mfit_object_path}/output/00{self.unseq[i]}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c_TAC.fits')
            #Run Molecfit if not already done
            if not self.filename[i].is_file() and self.redo_molec:
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
                    self.list_failed.append([self.night[i],self.unseq[i],e])
                    continue
                self.run_molecfit()
            try:
                self.wave[i], self.flux, self.flux_corr = self.load_FITS_table_data()
                logger.info('Successfully loaded Molecfit corrected spectrum')
            except Exception as e:
                logger.critical(f'Problem with the Molecfit corrected spectrum: {e}')
                self.list_failed.append([self.night[i],self.unseq[i],e])
                continue
            self.rv[i],self.err_rv[i] = self.calc_rad_vel()
            self.wave[i] = self.correct_vel_shift()
            self.model_flux[i] = self.load_model_SED()
            self.response[i] = self.get_response()
            # self.wl_norm[i] = self.get_norm_wl_array()
            # self.spline_fit[i] = self.fit_poly(self.poly_degree)
            self.spline_fit[i] = self.fit_spline()
            self.residuals[i] =  self.response[i] - self.spline_fit[i]
            # self.rms,self.bic = self.calc_goodness_of_fit(self.poly_degree)
            logger.info(f'Spline fit was successful.')
            # logger.debug(f'Residuals -  RMS : {self.rms}, BIC : {self.bic}')
            self.create_plots()
            self.save_response()
            logger.info(f'Response successfully saved to {self.rpath[i]}')
        if len(self.list_failed):
            failed = pd.DataFrame(data=self.list_failed)
            failed.to_csv(Path(f'/STER/karansinghd/PhD/ResponseCorrection/Failed/ListFailed_{self.night[i]}.txt'),index=False)
            # np.savetxt(Path(f'/STER/karansinghd/PhD/ResponseCorrection/Failed/ListFailed_{self.night[i]}.txt'),self.list_failed,fmt="%8d"+"%8d"+"%50s")
            logger.info('Saved list of failed corrections')

    def load_FITS_table_data(self):
        '''
        Function to load table data from the FITS file output by Molecfit.
        '''
        hdul = fits.open(self.filename[i])
        head = fits.getheader(self.filename[i])
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
        '''
        global i
        #load Balmer lines list
        # lines = np.loadtxt(self.Mask_path,usecols=(1),skiprows=1,delimiter =',',unpack=True)
        stype = spec_types[self.object[i]]
        lines = np.loadtxt(f'{self.LinesDir}/{fname[stype]}',delimiter =',',unpack=True,skiprows=1,usecols=[1])
        self.lines=lines
        #Create CCF
        specinterpol = si.interp1d(self.wave[i],self.flux,fill_value='extrapolate')
        velocity_array_full = np.array([])
        ccf_full = np.array([])
        #Look from -150 to +150 Km/s and get the CCF
        for k in np.arange(-150,150,0.25):
            linesred = np.add(lines,np.array(1000 * k/sc.c * lines))
            fluxred = specinterpol(linesred)
            velocity_array_full = np.append(velocity_array_full,float(k))
            ccf_full = np.append(ccf_full,fluxred.sum())

        #Get a small section around the minimum
        mini = ccf_full.argmin()
        if mini <= 20:
            limmin = 0
        else:
            limmin = mini-20

        if mini >= len(ccf_full) - 21:
            limmax = len(ccf_full) - 1
        else:
            limmax = mini + 20

        velocity_array = velocity_array_full[limmin:limmax]
        ccf = ccf_full[limmin:limmax]
        maxi = ccf_full.argmax()

        #Bisector analysis
        depth = ccf.max()-ccf.min()
        scale = depth/10
        bisector_y = np.linspace(ccf.min()+scale,ccf.max()-6*scale,20)

        minimumccf = ccf.argmin()
        ccfinterpol_blue = si.interp1d(ccf[0:minimumccf],velocity_array[0:minimumccf],fill_value='extrapolate')
        ccfinterpol_red = si.interp1d(ccf[minimumccf:],velocity_array[minimumccf:],fill_value='extrapolate')

        leftarray = ccfinterpol_blue(bisector_y)
        rightarray = ccfinterpol_red(bisector_y)

        velocityarray = (rightarray-leftarray)/2.+leftarray

        rv = velocityarray.mean()
        err_rv = velocityarray.std()
        logger.info(f'The topocentric velocity is {rv:.2f} km/s with a standard deviation of {err_rv:.2f} km/s')
        return rv,err_rv

    def correct_vel_shift(self):
        '''Function to correct for the RV shift before getting  the  response'''
        radvel = self.rv[i]*1000 #in m/s
        logger.info(f'Correcting for radial velocity of:{self.rv[i]:.2f} Km/s')
        return self.wave[i]*(1+radvel/sc.c)

    def load_model_SED(self):
        '''Function to load the model SED of the corresponding program 2 star'''
        logger.debug(f'Model Path: {self.InputFilesDir}/{ref_sed[self.object[i]]}')
        ref_wave, ref_flux = np.genfromtxt(fname=f'{self.InputFilesDir}/{ref_sed[self.object[i]]}', unpack=True)
        #interpolate the model at the wavelengths of the spectrum
        model_flux = si.interp1d(ref_wave, ref_flux, kind='linear',fill_value='extrapolate')(self.wave[i])
        logger.info('Successfully loaded model SED')
        return model_flux

    def get_response(self):
        '''Function to get the response and apply a median filter '''
        self.rough_response = self.flux_corr/self.model_flux[i]
        #2 median filters to deal with the edges
        R1_1 = mf(self.rough_response,1001)
        R1_2 = mf(self.rough_response,11)
        R1_1[:200] = R1_2[:200]
        R1_1[-40:] = R1_2[-40:]
        return R1_1

    def get_norm_wl_array(self):
        '''Function to normalize the array so that a polynomial fit of a high
        degree is possible. The array has to be sorted (it is by default)
        Wavelength array scaled between 1 and 3'''
        maximum = self.wave[i][-1]
        minimum = self.wave[i][0]
        factor = (maximum)/2
        return (((self.wave[i]-minimum)/factor)+1).astype(float)

    def fit_poly(self,d):
        '''Function to fit a polynomial to the median filtered response
        update (20201020): changed to spline fitting :) '''
        x=self.wl_norm[i]
        y=self.response[i]
        # d=self.poly_degree
        #Giving a higher weight to the edges improves the fit!
        weights = np.ones((len(y)))
        weights[2000:len(y)-2000] = 0.8
        c = poly.polyfit(x,y,d,w=weights)
        # c = poly.polyfit(x,y,d)
        polynomial = poly.polyval(x,c)
        logger.info(f'Fit the response with a polynomial of degree {self.poly_degree}')
        return polynomial

    def fit_spline(self):
        x=self.wave[i]
        y=self.response[i]

        left = min(x)+10
        right = max(x)-10

        #create knot array
        #[3781.,3852.,3913.,3952.,4000., 4056.,4112.,4168.,4224.,4280.
        added = np.array([3781,3852,3913,3952])
        violet = np.linspace(4000,4280,6)
        # violet2 = np.linspace(left,4280,10)
        blue = np.linspace(4350,5515,30)
        green = np.linspace(5545,6635,30)
        red = np.linspace(6675,right,20)

        knots = np.concatenate((added,violet,blue,green,red))
        self.knots=knots
        # knots2 = np.concatenate((violet2,blue,green,red))
        # print(knots)
        spl = us(x,y,t=knots)
        # spl2 = us(x,y,t=knots2)
        # f,(a,b,c) = pl.subplots(nrows=3,sharex=True,gridspec_kw={'height_ratios':[1,1,1]})
        # b.plot(x,y,'k')
        # [c.axvline(x=j,color='grey',alpha=0.8) for j in knots]
        # b.plot(x,spl(x),'b',label='fixed knots')
        # # b.plot(x,spl2(x),'r',ls='--',label='array knots')
        # b.legend()
        # # print(self.flux_corr)
        # # a.plot(self.wave[i],self.flux_corr)
        # c.plot(self.wave[i],self.flux_corr/spl(x),'k',lw=1.5,label='fixed knots')
        # # c.plot(self.wave[i],self.flux_corr/spl2(x),'r',ls='--',alpha=0.8,label='array knots',lw=0.8)
        # c.legend()
        # a.plot(self.wave[i],self.model_flux[i])
        # [b.axvline(x=j,color='grey',alpha=0.8) for j in knots]
        # [a.axvline(x=j,color='grey',alpha=0.8) for j in knots]
        # f.tight_layout()
        # c.set_ylabel('RC flux')
        # a.set_ylabel('Model')
        # b.set_ylabel('Response')
        # pl.show()
        # exit()
        return(spl(x))

    def calc_goodness_of_fit(self,d):
        '''Function to calc RMS and BIC value of the fit.
        '''
        logger.info('Calculating goodness of fit')

        y=self.residuals[i]
        N=len(y)

        rms = np.sqrt(np.sum(y**2)/N)
        BIC=rp.bic.bic(self.response[i],self.spline_fit[i],N)

        return rms,BIC


    def save_response(self):
        '''Function to save the response as a text file'''
        savefile=pd.DataFrame({'wavelength':self.wave[i],'spline':self.spline_fit[i],'response':self.response[i]})
        savefile.to_csv(self.rpath[i],header=True,index=False,sep='\t')

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
        with self.par_path[i].open(mode='w') as fout:
            with self.Mfit_dir.joinpath('template.par').open(mode="r") as f:
                lines=f.readlines()
                for line in lines:
                    if line.strip().startswith('user_workdir'):
                        line = f'user_workdir: {self.mfit_object_path}'
                        logger.debug(f'user_workdir updated: {self.mfit_object_path}')
                    elif line.strip().startswith('filename'):
                        line=f'filename: {self.mfit_object_path}/00{self.unseq[i]}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits'
                        logger.debug(f'filename updated: {self.mfit_object_path}/00{self.unseq[i]}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits')
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
                        line=f'output_name: {self.object[i]}_{self.unseq[i]}'
                        logger.debug(f'output name updated: {self.object[i]}_{self.unseq[i]}')
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
        activate_comm = sp.run(f'source /home/karansinghd/.bashrc; conda activate py27; /home/karansinghd/Applications/Molecfit/bin/molecfit {self.par_path[i]};/home/karansinghd/Applications/Molecfit/bin/calctrans {self.par_path[i]};conda deactivate',shell=True)
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
        ifn = self.HermesPath.joinpath(f'{self.night[i]}/reduced/00{self.unseq[i]}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits')
        ifn_old = self.HermesPath.joinpath(f'{self.night[i]}/reduced/{self.unseq[i]}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits')

        ifn_var = self.HermesPath.joinpath(f'{self.night[i]}/reduced/00{self.unseq[i]}_HRF_OBJ_ext_CosmicsRemoved_log_mergedVar_c.fits')
        ifn_var_old = self.HermesPath.joinpath(f'{self.night[i]}/reduced/{self.unseq[i]}_HRF_OBJ_ext_CosmicsRemoved_log_mergedVar_c.fits')
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
                logger.critical(f'Cannot find the spectrum (UNSEQ:{self.unseq[i]}, NIGHT:{self.night[i]}) in the HERMES database!')
        try:
            var = fits.getdata(ifn_var)
            self.flag_var=1
        except:
            try:
                var = fits.getdata(ifn_var_old)
                self.flag_var=1
            except:
                self.flag_var=0
                logger.critical(f'Cannot find the variance (UNSEQ:{self.unseq[i]}, NIGHT:{self.night[i]}) in the HERMES database!')
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
        ifn = self.HermesPath.joinpath(f'{self.night[i]}/reduced/00{self.unseq[i]}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits')
        ifn_old = self.HermesPath.joinpath(f'{self.night[i]}/reduced/{self.unseq[i]}_HRF_OBJ_ext_CosmicsRemoved_log_merged_c.fits')
        try:
            hdul = fits.open(ifn)
        except:
            try:
                hdul = fits.open(ifn_old)
            except:
                logger.critical(f'Cannot find the spectrum (UNSEQ:{self.unseq[i]}, NIGHT:{self.night[i]}) in the HERMES database!')
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
                # print('HRF' in df_good_models.filename[0])
                # print(df_good_models[['time_delta','night']])
                # df_good_models.drop(labels='time_delta',inplace=True)
                # print(df_good_models.columns)
                # exit()
                # df_good_models = df_good_models.iloc[[(df_good_models['night']-self.Night).abs().argmin()]]
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
        #Comment out next 3 lines if running the automatic loop!
        # save = input('Create and save plot ([y]/n)?\n')
        # if save.lower()=='n':
        #     return
        logger.info('Creating a plot of the fit')
        # mpl.rc('lines',lw=3)

        # self.flux_corr/self.model_flux[i]
        fig,(ax3,ax4,ax1,ax2)=pl.subplots(nrows=4, ncols=1, sharex=True, sharey=False, gridspec_kw={'height_ratios':[1,1,1,1]},figsize=(8,8))
        # fig,(ax1,ax2)=pl.subplots(nrows=2, ncols=1, sharex=True, sharey=False, gridspec_kw={'height_ratios':[3,1]},figsize=(8,8))
        ax3.plot(self.wave[i],self.flux_corr)
        ax4.plot(self.wave[i],self.model_flux[i])
        ax4.axvline(x = self.lines,color='cornflowerblue',alpha=0.9,ls='--')
        ax1.plot(self.wave[i],self.rough_response,'b',label='Rough response',alpha=0.5)
        ax1.plot(self.wave[i],self.response[i],'r',label='Median filtered response')
        ax1.plot(self.wave[i],self.spline_fit[i],'k--',label=f'Spline fit')
        ax2.plot(self.wave[i],self.residuals[i],'k')
        ax2.axhline(y=0,color='r')

        [ax1.axvline(x=j,color='grey',alpha=0.8) for j in self.knots]

        ax2.set_xlabel('Wavelength ($\AA$)')
        ax1.set_ylabel('ADU')
        # ax1.set_ylabel('ADU/(erg cm$^{-2}$ s$^{-1}$ $\AA^{-1}$)')
        ax2.set_ylabel('Residuals')

        ax1.set_title(f'NIGHT: {self.night[i]}, OBJECT:{self.object[i]}, UNSEQ:{self.unseq[i]}')
        ax1.legend(loc='upper right')
        fig.tight_layout()
        pl.savefig(f'/STER/karansinghd/PhD/Projects/P28_c/plots/{self.night[i]}_{self.object[i]}_{self.unseq[i]}.png',format='png')
        pl.show()
        pl.close()
        logger.info(f'Figure saved to /STER/karansinghd/PhD/Projects/P28_c/plots/{self.night[i]}_{self.object[i]}_{self.unseq[i]}.png')


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
    x = response(night=20101004,tolerance=None,overwrite=True)
    # correct_spectrum(20101220,325228,'HD36267')

    # HD36267, STDNIGHT: 20101220, STDUNUSEQ: 325226

    exit()
    # x=response(night=20190504,tolerance=None,overwrite=False)
    # x=response(20120917,tolerance=15,overwrite=True)
    df = pd.read_csv('melchiors_meta.csv',sep='|')
    df = df.drop_duplicates(subset=['night'])
    for _,row in df.iterrows():
        x = response(night=row['night'],tolerance=30,overwrite=True)
    # df.pop('Unnamed: 0')
    # df.pop('Unnamed: 31')
    # df.drop()
    print(df.head())
    exit()
    with Path('./response.log').open('w') as logfile:
        startdate=dt.date(year=2010,month=1,day=1)
        delta = dt.timedelta(days=1)
        enddate = dt.date(year=2011,month=12,day=31)
        print(startdate,enddate)
        while startdate <= enddate:
            try:
                x=response(night=int(startdate.strftime(format="%Y%m%d")),tolerance=None,overwrite=True)
                startdate += delta
            except Exception as e:
                logfile.write(f'{startdate.strftime(format="%Y%m%d")},{logger.critical(traceback.format_exc())}\n Error: {e}\n')
                print(e)
                startdate += delta
