################################
########### Imports ############
################################

import sys
import os
sys.path.append('../../python/')
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

#set up matplotlib
os.environ['MPLCONFIGDIR'] = '../mplstyles'
import matplotlib as mpl
from matplotlib import pyplot as plt
#got smarter about the mpl config: see mplstyles/ directory
plt.style.use('standard')

#fonts
# Set the font dictionaries (for plot title and axis titles)
title_font = {'fontname':'Arial', 'size':'16', 'color':'black', 'weight':'normal',
              'verticalalignment':'bottom'} # Bottom vertical alignment for more space
axis_font = {'fontname':'Arial', 'size':'32'}
legend_font = {'fontname':'Arial', 'size':'22'}

#fonts global settings
mpl.rc('font',family=legend_font['fontname'])

import numpy as np
import pandas as pd
from scipy.optimize import minimize
import emcee

#For some reason, this speeds things up!
os.environ["OMP_NUM_THREADS"] = "1"
################################
######Module Variables##########
################################
#see: https://docstore.mik.ua/orelly/other/python/0596001886_pythonian-chp-7-sect-1.html 
#for conventions on private variables (i.e. leading underscores)

refit=False
#fit parameters
_fitk=0.137
_fitq=1e-3

################################
########Fit Settings ###########
################################
#Construct a dictionary to store all the MCMC fit parameters and results
#These are all the settings a user will regularly change

        ########################## Data Settings ##########################
_mcmc_data={'g4_load_frac':1,
          'cap_load_frac':1,
          #'cap_sim_file':'/data/chocula/villaa/cascadeSimData/si28_R68_400k.pkl',
          #'cap_rcapture':0.161,
          'cap_sim_file':'../data/v3_400k.pkl',
          'cap_rcapture':0.218,
           ########################## Spectrum Settings ##########################
          'Emax': 2000, #[eVee]
          'Ebins': np.linspace(0,2500,251), #np.linspace(0,2000,201),
           'Efit_min':50, #[eVee]
           'Efit_max':2000, #1750, #[eVee]
           'spectrum_units':'reco-rate', #One of {'counts', 'reco-rate'}
           ########################## Yield Model Settings ##########################
           #'Y_model':'Lind',
           #'Y_labels': [r'k', r'$F_{NR}$'],
           #'Y_bounds': [(0.05,0.3),(0,30)],
           #'Y_model':'Chav',
           #'Y_labels': [r'k', r'$a^{-1}$', r'$F_{NR}$'],
           #'Y_bounds': [(0.05,0.3),(0,2e3),(0,30)],
           'Y_model':'Sor',
           'Y_labels': [r'k', r'q', r'$F_{NR}$'],
           'Y_bounds': [(0.05,0.3),(0,3e-2),(0,30)],
           #'Y_model':'AC',
           #'Y_labels': [r'k', r'$\xi$', r'$F_{NR}$'],
           #'Y_bounds': [(0.05,0.3),(0,2e3),(0,30)],
           #'Y_model':'pchip',
           #'Y_labels': [r'k', 'Er0', 'Er1', 'Er2', 'f1', r'$F_{NR}$'],
           #'Y_bounds': [(0.05,0.3),(0,1e-3),(0,1e-3),(0,1e-3),(0,1),(0,10)],
           #'Y_model':'Shexp',
           #'Y_labels': [r'k', 'Yshelf', 'Ec', 'dE', 'alpha', r'$F_{NR}$'],
           #'Y_bounds': [(0.05,0.3),(0,0.3),(0,1e3),(0,1e3),(0,100),(0,30)],
           #'Y_model':'Pol3',
           #'Y_labels': [r'p0', r'p1', r'p2', r'$F_{NR}$'],
           #'Y_bounds': [(-0.5,0.5),(-5e-4,5e-4),(-5e-7,5e-7),(0,30)],
           ########################## Sim Spectra Settings ##########################
           'ER_spec_model':'sim', #One of {'sim', 'flat') to selct from G4 simulation or flat
           #'ER_par_labels':[r'$scale_{G4}$'],
           'ER_par_labels':[r'$scale_{ER}$'],
           'ER_par_bounds':[(0,20)], #Unitless scaling factor
           #'ER_spec_model':'flat',
           #'ER_par_labels':[r'$R0_{ER}$'],
           #'ER_par_bounds':[(0,4e-2)], # Units are [Counts/sec/eVee bin] or [Counts/eVee bin] depending on spectrum_units
           #
           'NR_spec_model':'sim', #One of {'sim', 'flat', 'exp') to selct from G4 simulation, flat, or exponential
           #'NR_par_labels':[r'$scale_{G4}$'],
           'NR_par_labels':[r'$scale_{NR}$'],
           'NR_par_bounds':[(0,20)], #Unitless scaling factor
           #'NR_spec_model':'exp',
           #'NR_par_labels':[r'$R0_{NR}$',r'$E0_{NR}$'], #R0*exp(-E/E0) gives NR spectrum (post-yield)
           #'NR_par_bounds':[(0,0.1),(0,2e3)], # Units are [Counts/sec/eVee bin, eVee] or [Counts/eVee bin] depending on spectrum_units
           #
           'NG_spec_model':'sim', #Not going to implement anything other than sim for (n,gamma) yet
           'NG_par_labels':[r'$scale_{ng}$'],
           'NG_par_bounds':[(0,10)], #Unitless scaling factor
           ########################## Likelihood Settings ##########################
           'likelihood':'SNorm', #One of {'Pois', 'Norm', 'SNorm'} Only SNorm uses sigmas, others assume Pois stats
           ########################## Uncertainty Settings ##########################
           'doDetRes': True, #Include detector resolution effects
           'fpeak':1, #0.753 -- 1.0
           'doEffsyst':False, #Include systematics from cut efficiencies
           'doBurstLeaksyst':False, #Include burst cut leakage systematic
           ########################## MCMC Settings ##########################
           'nwalkers':128,
           'nstep':500000,
           'guesses':'Uniform', #Can either be uniform or shape (nwalkers, ndim),
           'moves':'DE8020',#'Default': StretchMove, 'DE8020': 80/20 DEMove/DESnookerMove
           'saveMCMC':True
          }


################################
###########Functions############
################################

def parse_options():
    _mcmc_data['labels']=_mcmc_data['Y_labels']+_mcmc_data['ER_par_labels']+_mcmc_data['NR_par_labels']+_mcmc_data['NG_par_labels']
    _mcmc_data['bounds']=_mcmc_data['Y_bounds']+_mcmc_data['ER_par_bounds']+_mcmc_data['NR_par_bounds']+_mcmc_data['NG_par_bounds']

    #Special case if ER and NR are both sim and we want to use the same G4 scaling factor for both:
    #if (mcmc_data['ER_spec_model']=='sim') and (mcmc_data['NR_spec_model']=='sim'):
    if (_mcmc_data['ER_par_labels']==[r'$scale_{G4}$']) and (_mcmc_data['NR_par_labels']==[r'$scale_{G4}$']):
      _mcmc_data['labels']=_mcmc_data['Y_labels']+_mcmc_data['NR_par_labels']+_mcmc_data['NG_par_labels']
      _mcmc_data['bounds']=_mcmc_data['Y_bounds']+_mcmc_data['NR_par_bounds']+_mcmc_data['NG_par_bounds']
    
    _mcmc_data['ndim']=len(_mcmc_data['labels'])

    return True

################################
######Execute on Load###########
################################
def printfunc: 
	print(_mcmc_data) #troubleshooting line
printfunc()
parse_options()
