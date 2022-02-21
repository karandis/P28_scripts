import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
from pathlib import Path
from scipy.constants import c
import scipy.interpolate as si

InputFilesDir = Path('/STER/karansinghd/PhD/ResponseCorrection/ModelSEDs')
LinesDir = Path('/STER/karansinghd/PhD/ResponseCorrection/VelocityCorrection/masks/')

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
# fname = {'A':'A_68Tau_depth0.95_blend20.list','B':'B_HR7512_depth0.98_blend20.list','F':"HD84937.txt"}
fname = {'A':'A_68Tau_depth0.95_blend20.list','B':'B_HR7512_depth0.98_blend20.list','F':"F_Procyon_depth0.95_blend20.list"}

def calc_rad_vel(wave,flux,stype,lines):

    #Create CCF
    specinterpol = si.interp1d(wave,flux,fill_value='extrapolate')
    velocity_array_full = np.array([])
    ccf_full = np.array([])
    #Look from -150 to +150 Km/s and get the CCF
    for k in np.arange(-150,150,0.25):
        linesred = np.add(lines,np.array(1000 * k/c * lines))
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
    return rv,err_rv



starnames = ['HD152614','HD36267','HD118098','HD14055','HD87887','HD46300', 'HD184006','HD149212','HD147449','GSC4293-0432','HD214994','HD42818','HD56169','HD185395','HD206826','HD220657','HD84937']

for star in starnames:
    star = 'HD56169'
    f,a = pl.subplots()

    ref_wave, ref_flux = np.genfromtxt(fname=f'{InputFilesDir}/{ref_sed[star]}', unpack=True)
    #introduce RV shift
    # radvel = -8000 #in m/s
    # ref_wave*=(1+radvel/c)
    #####
    stype = spec_types[star]
    # if stype =='F' or stype == 'A':
        # continue
    '''
    lines = np.loadtxt(f'{LinesDir}/{fname[stype]}',delimiter =',',unpack=True,skiprows=1,usecols=[1])

    df = pd.read_csv(f'{LinesDir}/{fname[stype]}',header=0,sep=',')
    lines = df.loc[df.loggf >=0.5].wavelength.to_numpy()

    lines = np.loadtxt('linelists/HD220657.txt')
    '''

    try:
        lines = np.loadtxt(f'{LinesDir}/{fname[stype]}',unpack=True)
    except:
        lines = np.loadtxt(f'{LinesDir}/{fname[stype]}',delimiter =',',unpack=True,skiprows=1,usecols=[1])

    #find shifted RV
    rv,err = calc_rad_vel(ref_wave,ref_flux,stype,lines)
    print(rv,err)
    # exit()
    a.plot(ref_wave,ref_flux,'k')
    [a.axvline(x=line,color='k',ls='--') for line in lines]
    a.set_title(f'{star}')
    a.set_xlim(min(lines)-5,max(lines)+5)
    pl.show()
    exit()
