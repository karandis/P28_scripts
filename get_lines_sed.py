import numpy as np
import matplotlib.pyplot as pl

from pathlib import Path

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
fname = {'A':'A_68Tau_depth0.95_blend20.list','B':'B_HR7512_depth0.98_blend20.list','F':"F_Procyon_depth0.95_blend20.list"}

starnames = ['HD152614','HD36267','HD118098','HD14055','HD87887','HD46300', 'HD184006','HD149212','HD147449','GSC4293-0432','HD214994','HD42818','HD56169','HD185395','HD206826','HD220657','HD84937']

for star in starnames:
    f,a = pl.subplots()

    ref_wave, ref_flux = np.genfromtxt(fname=f'{InputFilesDir}/{ref_sed[star]}', unpack=True)
    stype = spec_types[star]
    if stype =='F' or stype == 'A':
        continue
    lines = np.loadtxt(f'{LinesDir}/{fname[stype]}',delimiter =',',unpack=True,skiprows=1,usecols=[1])
    a.plot(ref_wave,ref_flux,'k')
    [a.axvline(x=line,color='k',ls='--') for line in lines]
    a.set_title(f'{star}')
    a.set_xlim(min(lines)-5,max(lines)+5)
    pl.show()
    exit()
