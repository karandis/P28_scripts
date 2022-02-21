directory of scripts for the MELCHIORS project (P28)

correct_spectrum*.py - scripts to correct the spectra for the HERMES instrumental response

get_response*.py - script to calculate the instrumental response, usually called by correct_spectrum.py. Can also be used as a standalone script.

get_response_gienah.py, get_response_old.py, get_response.py - poor version control remnants

README

I. Overview of files

'compare_reductions.py'
- script to compare the processed P28 data from April 2021 vs Jan 2022. Creates 'modified.csv'

'modified.csv'
- List of P28 stars which whose 'stdinfo.csv' entries were modified with the new reduction in Jan 2022. Created with 'compare_reductions.py'

'correct_spectrum.py'
- script to correct spectra for the instrumental response. Relies on the response already being derived in '/STER/karansinghd/PhD/ResponseCorrection/responses_c_2022/'. If no response exists for the required night, the script calls 'get_response.py'.

'diagnostic_plots.py'
- Script to create diagnostic plots starting from '_c.fits' until the final '_cr.fits'.

'fix_fits.py'
- Script to fix the headers for a subset of '/STER/karansinghd/P28_cr_2022/stdinfo.csv'. Something went wrong when Molecfit saved header["UTC"], and so this had to be redone.

"FunkyBlueData.csv"
- File containing list of FunkyBlue spectra from '/STER/pierre/hermes/p28/funkyblue', with auxillary information.

'list_bad_spec.py'
- Script to create 'FunkyBlueData.csv' and some plots to check if there was some correlation with airmass, P2 star, etc.

'get_response.py'
- Main script to get the instrumental response of a night, given a particular tolerance. The response is saved in '/STER/karansinghd/PhD/ResponseCorrection/responses_c_2022/'.

'Make_plots_molecfit.py'
'Make_plots.py'
- Scripts to create plots used in the paper.

'py27.yaml'
- The yaml file for my python 2.7 environment with packages that will allow the user to run the standalone version of Molecfit (see below)

'plotting.py'
- List of matplotlib presets. I have added this script to my PYTHONPATH but you can also keep this in the same directory. Delete or modify to your use!

II. Prerequisites

- The software MOLECFIT (https://www.eso.org/sci/software/pipelines/skytools/molecfit).

  The standalone versions are available via the FTP on the provided link (https://ftp.eso.org/pub/dfs/pipelines/skytools/molecfit/). However, they are no longer supported.

  The new package is downloadable via the pipeline (https://www.eso.org/sci/software/pipelines/#pipelines_table). NOTE: I have not tried this, I only used the standalone version since it was the only option at the time!

- A python 2.7 environment that can run Molecfit. You can install the environment using 'conda env create -f py27.yml'. Note, this might be different with the ESO reflex pipeline!

- A python 3.6 or later environment to run the other scripts. Basic modules required are 'astropy', 'pandas', 'numpy', 'matplotlib' and 'pathlib'. I also include the 'plotting.py' matplotlib presets that I have added to my PYTHONPATH, so I can import it without bothering about aesthetics.

III. Running the HERMES instrumental response correction

'get_response.py' is the main script needed to obtain the response for a night

E.g. Snippet to get the response for the night of Jan 20, 2022 with a tolerance of 60 days:

'''
import get_response as gr
resp = gr.response(night=20210120, tolerance=60, overwrite = True)
'''

The reason it is in the form of a class is that 'resp' has the wavelength, P2 flux, Molecfit corrected spectrum, model SED, response and the spline fit to the response. Some of these can be plotted as follows:

'''
import matplotlib.pyplot as pl
f,(a1,a2) = pl.subplots(nrows=2,ncols=1,sharex=True,gridspec_kw={'height_ratios':[4,1]},ffigsize=(8,6))
a1.plot(resp.wave,resp.response,'r')
a1.plot(resp.wave,resp.spline_fit,'k')
a2.plot(resp.wave,resp.residuals,'k')
a2.axhline(y=0,color='r')
pl.show()
'''

The list of available class variables is commented in 'get_response.py'

IV. Correcting a spectrum for the instrumental response

'correct_spectrum.py'

Scripts to correct HERMES spectra for the instrumental response. Can be used as follows:

'''
import correct_spectrum as cs
f = cs.correct_spectrum(night=20100930,unseq=307479,object='V__AG_Psc')
'''

The script will correct this particular spectrum (with and without Molecfit) and save the output to a preset directory that is hardcoded '/STER/karansinghd/PhD/Projects/P28_c', and then copy it to '/STER/karansinghd/P28_cr_2022'.


Both the response determination and the spectrum correction can be integrated into the HERMES pipeline with minimal effort.
