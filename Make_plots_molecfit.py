'''
Script to make plots for the paper

Creates a script showing the Molecfit correction in detail
'''

import numpy as np
import pandas as pd
import plotting
from pathlib import Path
import matplotlib.pyplot as pl
from astropy.io import fits
plotting.set_mode('paper')
'''
To check:
307479,V__AG_Psc
307530,307543,307614,363184
'''

s='/STER/karansinghd/PhD/Projects/P28_c/V__AG_Psc/output/00307479_HRF_OBJ_ext_CosmicsRemoved_log_merged_c_TAC.fits'

hdu = fits.open(s)
table = hdu[1].data
w = table.field(0)
f = table.field(2)
f_c = table.field(5)
print(table)
ix = (w>=6200) & (w<9000)
w=w[ix]
f=f[ix]
f_c=f_c[ix]
f/=1000
f_c/=1000

fig,ax = pl.subplots()
ax.plot(w,f,'r',label='Original spectrum')
ax.plot(w,f_c,'k',label='Molecfit corrected')

ax.set_xlim([6200,9002])
ax.set_xlabel('Wavelength ($\AA$)')
ax.set_ylabel('Flux (ADU)')
fig.tight_layout()
ax.legend()
# fig.savefig('/Users/karansinghd/PhD/Projects/P28/Molecfit_correction.pdf',dpi=300,format='pdf',bbox_inches='tight')
pl.show()
