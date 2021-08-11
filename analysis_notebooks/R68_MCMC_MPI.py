#This is a script to run MCMC fits using MPI
#It is roughly a script version of R68_MCMC(_faster).ipynb

#Run using command like:
#mpiexec --hostfile /home/mastx027/MPI/mpi_hostfile -n 24 python R68_MCMC_MPI.py > out.txt 2>&1 &
#Make sure you're in the nr_fano conda environment

################################################################################
#Setup
################################################################################

#we may need some code in the ../python directory and/or matplotlib styles
import sys
import os
sys.path.append('../python/')
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
from constants import *

#For some reason, this speeds things up!
os.environ["OMP_NUM_THREADS"] = "1"



################################################################################
#Fit Settings
################################################################################
#Construct a dictionary to store all the MCMC fit parameters and results
#These are all the settings a user will regularly change

        ########################## Data Settings ##########################
mcmc_data={'g4_load_frac':0.1,
          'cap_load_frac':0.1,
          #'cap_sim_file':'/data/chocula/villaa/cascadeSimData/si28_R68_400k.pkl',
          #'cap_rcapture':0.161,
          'cap_sim_file':'data/v3_400k.pkl',
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
           #'Y_model':'Sor',
           #'Y_labels': [r'k', r'q', r'$F_{NR}$'],
           #'Y_bounds': [(0.05,0.3),(0,3e-2),(0,30)],
           'Y_model':'Sor',
           'Y_labels': [r'k', r'$\xi$', r'$F_{NR}$'],
           'Y_bounds': [(0.05,0.3),(0,2e3),(0,30)],
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
           'fpeak':0.753, #0.753 -- 1.0
           'doEffsyst':True, #Include systematics from cut efficiencies
           'doBurstLeaksyst':True, #Include burst cut leakage systematic
           ########################## MCMC Settings ##########################
           'nwalkers':128,
           'nstep':500000,
           'guesses':'Uniform', #Can either be uniform or shape (nwalkers, ndim),
           'moves':'DE8020',#'Default': StretchMove, 'DE8020': 80/20 DEMove/DESnookerMove
           'saveMCMC':True
          }

################################################################################
#Parse Options
################################################################################
mcmc_data['labels']=mcmc_data['Y_labels']+mcmc_data['ER_par_labels']+mcmc_data['NR_par_labels']+mcmc_data['NG_par_labels']
mcmc_data['bounds']=mcmc_data['Y_bounds']+mcmc_data['ER_par_bounds']+mcmc_data['NR_par_bounds']+mcmc_data['NG_par_bounds']

#Special case if ER and NR are both sim and we want to use the same G4 scaling factor for both:
#if (mcmc_data['ER_spec_model']=='sim') and (mcmc_data['NR_spec_model']=='sim'):
if (mcmc_data['ER_par_labels']==[r'$scale_{G4}$']) and (mcmc_data['NR_par_labels']==[r'$scale_{G4}$']):
    mcmc_data['labels']=mcmc_data['Y_labels']+mcmc_data['NR_par_labels']+mcmc_data['NG_par_labels']
    mcmc_data['bounds']=mcmc_data['Y_bounds']+mcmc_data['NR_par_bounds']+mcmc_data['NG_par_bounds']
    
mcmc_data['ndim']=len(mcmc_data['labels'])
################################################################################
#Load Data
################################################################################
import R68_load as r68
import R68_spec_tools as spec
#Import likelihood functions
from likelihoods import *
import R68_efficiencies as eff

#Set eVee energy binning
Emax=mcmc_data['Emax']
Ebins=mcmc_data['Ebins']
Ebins_ctr=(Ebins[:-1]+Ebins[1:])/2

#Set fit range
E_lim_min=mcmc_data['Efit_min'] #eVee
E_lim_max=mcmc_data['Efit_max'] #eVee
spec_bounds=(np.digitize(E_lim_min,Ebins)-1,np.digitize(E_lim_max,Ebins)-1)
mcmc_data['spec_bounds']=spec_bounds

#Measured spectra background subtraction
meas=r68.load_measured(verbose=True)
tlive_PuBe=meas['PuBe']['tlive']

#Note that N_meas may be either in counts or reconstructed rate, depending on mcmc_data['spectrum_units']
#uncertainty is (high,low)
N_meas,dN_meas=spec.doBkgSub(meas, Ebins, mcmc_data['Efit_min'], mcmc_data['Efit_max'],\
                             doEffsyst=mcmc_data['doEffsyst'], doBurstLeaksyst=mcmc_data['doBurstLeaksyst'],\
                             output=mcmc_data['spectrum_units'])

if mcmc_data['likelihood']=='SNorm':
    #Precalculate split normal likelihood params if we're going to need it
    SNpars=getSNparsArray(N_meas[slice(*spec_bounds)],dN_meas[0][slice(*spec_bounds)],dN_meas[1][slice(*spec_bounds)])
    SNpars=SNpars.T

#Load g4 Simulations
#No need for this if we're using analytical models for both
if (mcmc_data['ER_spec_model']=='sim') or (mcmc_data['NR_spec_model']=='sim'):
    g4=r68.load_G4(load_frac=mcmc_data['g4_load_frac'])

    #Trim events that won't figure into the analysis range
    #Trimmed sim data
    if (mcmc_data['ER_spec_model']=='sim'):
        Eee_er=np.sum(g4['ER']['E'],axis=1)
        Evec_er_cut=(Eee_er>10) & (Eee_er<3e3)
        Evec_er=g4['ER']['E'][Evec_er_cut]

    if (mcmc_data['NR_spec_model']=='sim'):
        Eee_nr=np.sum(g4['NR']['E'],axis=1)
        Evec_nr_cut=(Eee_nr>10) & (Eee_nr<30e3)
        Evec_nr=g4['NR']['E'][Evec_nr_cut]


#Load Capture Simulations
cap=r68.load_simcap(file=mcmc_data['cap_sim_file'], rcapture=mcmc_data['cap_rcapture'], load_frac=mcmc_data['cap_load_frac'])

#Pre-calculate Simulated ER spectrum
#This is independent of other params, so we can do this once and reuse it
if mcmc_data['ER_spec_model']=='sim':
    if mcmc_data['spectrum_units']=='counts':
        N_er = spec.buildAvgSimSpectrum_ee(Ebins=Ebins, Evec=Evec_er, Yield=1.0, F=F, scale=1,\
                                           doDetRes=mcmc_data['doDetRes'], fpeak=mcmc_data['fpeak'],\
                                           doEffs=True)    
    elif mcmc_data['spectrum_units']=='reco-rate':
        N_er = spec.buildAvgSimSpectrum_ee(Ebins=Ebins, Evec=Evec_er, Yield=1.0, F=F, scale=1,\
                                           doDetRes=mcmc_data['doDetRes'], fpeak=mcmc_data['fpeak'],\
                                           doEffs=False)
elif mcmc_data['ER_spec_model']=='flat':
    N_er = np.ones_like(Ebins_ctr)
else:
    print('ER_spec_model: ',mcmc_data['ER_spec_model'],' not implemented yet.')
        
#Import yield models
import R68_yield as Yield
#Initialize Yield model
Y=Yield.Yield('Lind',[0.15])
model=mcmc_data['Y_model']
Y=Yield.Yield(model,np.zeros(Y.model_npar[model]))

################################################################################
#Calculate Log probability, log(likelihood*prior)
################################################################################
#theta: array of fit parameters (yield_par0, yield_par1, ...,  F_NR, scale_g4, scale_ng, ...)
#    These are in the same order as mcmc_data['labels']
#theta_bounds: paramter bounds, shape should be len(theta)x2
#spec_bounds: range of bin numbers in spectrum to consider. The analysis range is [bin_low,bin_high)
#likelihood: Likelihood function, either 'Pois' or 'Norm'

def calc_log_prob(theta=[0.2, 1, 1, 1], theta_bounds=((0,1),(0,10),(0,10),(0,10)), spec_bounds=(5,101),
                  likelihood='Pois'):

    
    
    #Access the global data
    #These must be already defined!!!
    global mcmc_data, eff, N_meas, N_er, tlive_PuBe, Evec_nr, cap, Y
    
    ############
    #Parse fit params
    nYpar=Y.npars

    Y.set_pars(theta[:nYpar])
    F_NR=theta[nYpar]
    
    #Calculate the (log)prior first since we may not need to calculate the likelihood
    lp=lp_flat(theta, theta_bounds)
    if not np.isfinite(lp):
        return -np.inf
    
    #Special check for pchip yield
    #Will fail to solve if the parameters are not ordered.
    #This effectively includes some parameter ordering in their priors
    if not Y.solve():
        return -np.inf
        
    #Grab correct ER scaling factor from theta
    if mcmc_data['ER_spec_model']=='sim':
        if '$scale_{G4}$' in mcmc_data['labels']:
            scale_er=theta[mcmc_data['labels'].index('$scale_{G4}$')]/g4['ER']['tlive']
        else:
            scale_er=theta[mcmc_data['labels'].index('$scale_{ER}$')]/g4['ER']['tlive']
    elif mcmc_data['ER_spec_model']=='flat':
        scale_er=theta[mcmc_data['labels'].index(r'$R0_{ER}$')]
        
    ##########
    #Build the spectra
    if mcmc_data['spectrum_units']=='counts':
        #This includes detector resolution, triggering, and cut efficiency effects
        #Does NOT include livetime or write efficiency

        #NR
        if mcmc_data['NR_spec_model']=='sim':
            N_nr=spec.buildAvgSimSpectrum_ee(Ebins=Ebins, Evec=Evec_nr, Yield=Y, F=F_NR, scale=1,\
                                             doDetRes=mcmc_data['doDetRes'], fpeak=mcmc_data['fpeak'],\
                                             doEffs=True)
            if '$scale_{G4}$' in mcmc_data['labels']:
                scale_nr=theta[mcmc_data['labels'].index('$scale_{G4}$')]/g4['NR']['tlive']
            else:
                scale_nr=theta[mcmc_data['labels'].index('$scale_{NR}$')]/g4['NR']['tlive']
        elif mcmc_data['NR_spec_model']=='exp':
            N_nr=np.exp(-Ebins_ctr/theta[mcmc_data['labels'].index(r'$E0_{NR}$')])
            scale_nr=theta[mcmc_data['labels'].index(r'$R0_{NR}$')]
            
        #(n,gamma)
        N_ng=spec.buildAvgSimSpectrum_ee_composite(Ebins=Ebins, Evec=cap['E'], dEvec=cap['dE'],\
                                                   Yield=Y, F=F_NR, scale=1, doDetRes=mcmc_data['doDetRes'],\
                                                   fpeak=mcmc_data['fpeak'], doEffs=True)
        scale_ng=theta[mcmc_data['labels'].index('$scale_{ng}$')]/cap['tlive']

        #Adjust for livetime and write efficiency
        N_pred = (N_nr*scale_nr + 
                  N_er*scale_er + 
                  N_ng*scale_ng)*tlive_PuBe*eff.eff_write
    
    elif mcmc_data['spectrum_units']=='reco-rate':
        #Don't apply any efficiency effects to simulated spectrum
        #NR
        if mcmc_data['NR_spec_model']=='sim':
            N_nr=spec.buildAvgSimSpectrum_ee(Ebins=Ebins, Evec=Evec_nr, Yield=Y, F=F_NR, scale=1,\
                                             doDetRes=mcmc_data['doDetRes'], fpeak=mcmc_data['fpeak'],\
                                             doEffs=False)
            if '$scale_{G4}$' in mcmc_data['labels']:
                scale_nr=theta[mcmc_data['labels'].index('$scale_{G4}$')]/g4['NR']['tlive']
            else:
                scale_nr=theta[mcmc_data['labels'].index('$scale_{NR}$')]/g4['NR']['tlive']
        elif mcmc_data['NR_spec_model']=='exp':
            N_nr=np.exp(-Ebins_ctr/theta[mcmc_data['labels'].index(r'$E0_{NR}$')])
            scale_nr=theta[mcmc_data['labels'].index(r'$R0_{NR}$')]
            
        #(n,gamma)
        N_ng=spec.buildAvgSimSpectrum_ee_composite(Ebins=Ebins, Evec=cap['E'], dEvec=cap['dE'],\
                                                   Yield=Y, F=F_NR, scale=1, doDetRes=mcmc_data['doDetRes'],\
                                                   fpeak=mcmc_data['fpeak'], doEffs=False)
        scale_ng=theta[mcmc_data['labels'].index('$scale_{ng}$')]/cap['tlive']
        
        #Calculate rate (though we'll still call them N_* just to be confusing)
        N_pred = (N_nr*scale_nr + 
                  N_er*scale_er + 
                  N_ng*scale_ng)

    ##########
    #Calculate the log probability = log prior + log likelihood
    ll=None
    
    if likelihood=='Norm':
        ll = ll_norm(N_meas[slice(*spec_bounds)],N_pred[slice(*spec_bounds)])
    elif likelihood=='Pois':
        ll = ll_pois(N_meas[slice(*spec_bounds)],N_pred[slice(*spec_bounds)])
    elif likelihood=='SNorm':
        ll = ll_SNorm(N_pred[slice(*spec_bounds)],*SNpars)
    else:
        print('Error: Bad likelihood')
        return None
    
    if not np.isfinite(ll):
        return -np.inf
    
    return lp + ll

################################################################################
#Fit Setup and Helper Function
################################################################################

import emcee
from multiprocessing import Pool

def Fit_helper(theta):                 
    return calc_log_prob(theta=theta, theta_bounds=tuple(mcmc_data['bounds']),
                         spec_bounds=mcmc_data['spec_bounds'], likelihood=mcmc_data['likelihood'])

################################################################################
#MCMC Fitting loop
################################################################################

import time
from schwimmbad import MPIPool

nwalkers=mcmc_data['nwalkers']
ndim=mcmc_data['ndim']
nstep=mcmc_data['nstep']

#Sample priors uniformly
if mcmc_data['guesses']=='Uniform':
    bounds_low=np.array(mcmc_data['bounds'])[:,0]
    bounds_hi=np.array(mcmc_data['bounds'])[:,1]
    mcmc_data['guesses']=(bounds_hi-bounds_low)*np.random.random_sample((nwalkers, ndim))+bounds_low
else:
    print('guesses = ',mcmc_data['guesses'],' not yet implemented')

with MPIPool() as pool:
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
    
    if mcmc_data['moves']=='Default':
        sampler = emcee.EnsembleSampler(nwalkers, ndim, Fit_helper, pool=pool)
    elif mcmc_data['moves']=='DE8020':
        sampler = emcee.EnsembleSampler(nwalkers, ndim, Fit_helper, pool=pool,\
                                        moves=[(emcee.moves.DEMove(), 0.8),(emcee.moves.DESnookerMove(), 0.2)])
    else:
        print('Error: mcmc move type not specified')
        
    start = time.time()
    sampler.run_mcmc(mcmc_data['guesses'], nstep, progress=True);
    end = time.time()
    mcmc_data['t_run']=end - start
    
    print('Total run time was {0:.2f} hours.'.format((end - start)/3600))
    
mcmc_data['sampler']=sampler
################################################################################
#Save this work
################################################################################

import pickle as pkl
import os
import misc

if mcmc_data['saveMCMC']:
    ifile = 0
    #fname='data/mcmc_{0}_{1}walk_{2}step_{3}_v{4}.pkl'.format(model,nwalkers,nstep,mcmc_data['likelihood'],ifile+1)
    fname='data/mcmc_{0}_{1}walk_{2}step_{3}_v{4}.pkl'.format(model,nwalkers,misc.human_format(nstep),
                                                              mcmc_data['likelihood'],ifile+1)
    while os.path.exists(fname):
        ifile += 1
        fname='data/mcmc_{0}_{1}walk_{2}step_{3}_v{4}.pkl'.format(model,nwalkers,misc.human_format(nstep),
                                                                  mcmc_data['likelihood'],ifile+1)
        
    print(fname)
    with open(fname, 'wb') as saveFile:
        pkl.dump(mcmc_data,saveFile)