#This is a script to run MCMC fits using MPI
#It is roughly a script version of R68_MCMC(_faster).ipynb

#Run using command like:
#mpiexec --hostfile /home/mast/MPI/mpi_hostfile -n 24 python R68_MCMC_MPI.py > out.txt 2>&1 &
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
mcmc_data={'g4_load_frac':0.1,
          'cap_load_frac':0.1,
          'cap_sim_file':'/data/chocula/villaa/cascadeSimData/si28_R68_400k.pkl',
          'cap_rcapture':0.161,
          'Emax': 2000, #[eVee]
          'Ebins': np.linspace(0,2000,201),
           'Efit_min':90, #[eVee]
           'Efit_max':1750, #[eVee]
           #'Ymodel':'Lind',
           #'labels': [r'k', r'$F_{NR}$', r'$scale_{G4}$', r'$scale_{ng}$'],
           #'theta_bounds': ((0.05,0.3),(0,30),(0.1,10),(0.1,10)),
           #'Ymodel':'Chav',
           #'labels': [r'k', r'$a^{-1}$', r'$F_{NR}$', r'$scale_{G4}$', r'$scale_{ng}$'],
           #'theta_bounds': ((0.05,0.3),(0,1e3),(0,30),(0.1,10),(0.1,10)),
           'Ymodel':'Sor',
           'labels': [r'k', r'q', r'$F_{NR}$', r'$scale_{G4}$', r'$scale_{ng}$'],
           'theta_bounds': ((0.05,0.3),(0,3e-2),(0,30),(0.1,10),(0.1,10)),
           #'Ymodel':'AC',
           #'labels': [r'k', r'$\xi$', r'$F_{NR}$', r'$scale_{G4}$', r'$scale_{ng}$'],
           #'theta_bounds': ((0.05,0.3),(0,2e3),(0,30),(0.1,10),(0.1,10)),
           #'Ymodel':'pchip',
           #'labels': [r'k', 'Er0', 'Er1', 'Er2', 'f1', r'$F_{NR}$', r'$scale_{G4}$', r'$scale_{ng}$'],
           #'theta_bounds': ((0.05,0.3),(0,1e-3),(0,1e-3),(0,1e-3),(0,1),(0,10),(0.1,10),(0.1,10)),
           #'Ymodel':'Shexp',
           #'labels': [r'k', 'Yshelf', 'Ec', 'dE', 'alpha', r'$F_{NR}$', r'$scale_{G4}$', r'$scale_{ng}$'],
           #'theta_bounds': ((0.05,0.3),(0,0.3),(0,1e3),(0,1e3),(0,100),(0,30),(0.1,10),(0.1,10)),
           'likelihood':'Pois',
           'doDetRes': True,
           'fpeak':1.0,
           'nwalkers':128,
           #'ndim':4,
           'ndim':5,
           #'ndim':8,
           'nstep':5000,
           'guesses':'Uniform', #Can either be uniform or shape (nwalkers, ndim),
           'moves':'DE8020',#'Default': StretchMove, 'DE8020': 80/20 DEMove/DESnookerMove
           'saveMCMC':True
          }

################################################################################
#Load Data
################################################################################

#Load the datasets
import R68_load as r68

meas=r68.load_measured()
g4=r68.load_G4(load_frac=mcmc_data['g4_load_frac'])
cap=r68.load_simcap(file=mcmc_data['cap_sim_file'], rcapture=mcmc_data['cap_rcapture'], load_frac=mcmc_data['cap_load_frac'])

#Import yield models
import R68_yield as Yield
import R68_spec_tools as spec

#Initialize Yield model
Y=Yield.Yield('Lind',[0.15])

model=mcmc_data['Ymodel']
Y=Yield.Yield(model,np.zeros(Y.model_npar[model]))

#Set eVee energy binning
Emax=mcmc_data['Emax']
Ebins=mcmc_data['Ebins']

#Measured spectra
N_meas_PuBe,_ = np.histogram(meas['PuBe']['E'],bins=Ebins)
N_meas_Bkg,_ = np.histogram(meas['Bkg']['E'],bins=Ebins)

tlive_PuBe = meas['PuBe']['tlive']
tlive_Bkg = meas['Bkg']['tlive']
#We'll scale everything to the PuBe live time and work with counts, not rate, to get the Poisson stats right

N_meas_Bkg_scaled = N_meas_Bkg * tlive_PuBe/tlive_Bkg
#Estimate of counts due to PuBe
N_meas = N_meas_PuBe - N_meas_Bkg_scaled

#g4 Simulations
#Trim events that won't figure into the analysis range
#Trimmed sim data
Eee_er=np.sum(g4['ER']['E'],axis=1)
Evec_er_cut=(Eee_er>10) & (Eee_er<3e3)
Evec_er=g4['ER']['E'][Evec_er_cut]

Eee_nr=np.sum(g4['NR']['E'],axis=1)
Evec_nr_cut=(Eee_nr>10) & (Eee_nr<30e3)
Evec_nr=g4['NR']['E'][Evec_nr_cut]

#Calculate Simulated ER spectrum
#This is independent of other params, so we can do this once and reuse it
N_er = spec.buildAvgSimSpectrum_ee(Ebins=Ebins, Evec=Evec_er, Yield=1.0, F=F, scale=1,\
                                   doDetRes=mcmc_data['doDetRes'], fpeak=mcmc_data['fpeak'])    

################################################################################
#Define likelihood functions
################################################################################
from scipy.special import factorial, gamma, loggamma

#Poisson likelihood of measuring k given expected mean of lambda
def pois_likelihood(k, lamb):
    return (lamb**k)*np.exp(-lamb)/gamma(k+1.)

#Poisson log-likelihood
#k: observed counts
#lamb: expected (model) counts
def ll_pois(k, lamb):   
    if np.sum(lamb<=0):
        return -np.inf
    
    return np.sum(k*np.log(lamb) - lamb - loggamma(k+1.))

#Normal log-likelihood, limit of Poisson for large lambda
#k: observed counts
#lamb: expected (model) counts
def ll_norm(k,lamb):
    if np.sum(lamb<=0):
        return -np.inf
    
    return np.sum(-0.5*np.log(2*np.pi*lamb) - (k-lamb)**2/(2*lamb))

#Log of flat prior functions
#theta: array of parameter values
#bounds: array of parameter bounds. shape should be len(theta)x2
def lp_flat(theta, bounds):
    #for itheta,ibounds in zip(theta,bounds):
    #    if not (ibounds[0] < itheta < ibounds[1]):
    #        return -np.inf
        
    #return 0.0
    
    if (np.array(bounds)[:,0]<theta).all() and (theta<np.array(bounds)[:,1]).all():
        return 0.0
    return -np.inf

#Log of normal prior distribution
#theta: parameter value(s)
#mu: parameter prior distribution mean(s)
#sigma: paramter prior distribution sigma(s)
def lp_norm(theta, mu, sigma):
    return np.sum(-0.5*((theta-mu)/sigma)**2 - np.log(sigma)-0.5*np.log(2*np.pi))


################################################################################
#Calculate Log probability, log(likelihood*prior)
################################################################################
#theta: array of fit parameters (yield_par0, yield_par1, ...,  F_NR, scale_g4, scale_ng, ...)
#theta_bounds: paramter bounds, shape should be len(theta)x2
#spec_bounds: range of bin numbers in spectrum to consider. The analysis range is [bin_low,bin_high)
#likelihood: Likelihood function, either 'Pois' or 'Norm'

def calc_log_prob(theta=[0.2, 1, 1, 1], theta_bounds=((0,1),(0,10),(0,10),(0,10)), spec_bounds=(5,101),
                  likelihood='Pois'):

    #Access the global data
    #These must be already defined!!!
    global N_meas, N_er, tlive_PuBe, Evec_nr, cap, Y
    
    ############
    #Set some local variables
    nYpar=Y.npars

    Y.set_pars(theta[:nYpar])
    F_NR=theta[nYpar]
    scale_g4=theta[nYpar+1]
    scale_ng=theta[nYpar+2]
    
    
    #Calculate the (log)prior first since we may not need to calculate the likelihood
    lp=lp_flat(theta, theta_bounds)
    if not np.isfinite(lp):
        return -np.inf
    
    #Special check for pchip yield
    #Will fail to solve if the parameters are not ordered.
    #This effectively includes some parameter ordering in their priors
    if not Y.solve():
        return -np.inf
        
    
    ##########
    #Build the spectra

    #NR
    N_nr=spec.buildAvgSimSpectrum_ee(Ebins=Ebins, Evec=Evec_nr, Yield=Y, F=F_NR, scale=1, 
                                               doDetRes=mcmc_data['doDetRes'], fpeak=mcmc_data['fpeak'])
    #(n,gamma)
    N_ng=spec.buildAvgSimSpectrum_ee_composite(Ebins=Ebins, Evec=cap['E'], dEvec=cap['dE'], Yield=Y, F=F_NR, scale=1, 
                                               doDetRes=mcmc_data['doDetRes'], fpeak=mcmc_data['fpeak'])
    

    #Total counts for PuBe live time
    #Uncertainty will be ~sqrt(N)
    N_pred = (N_nr*scale_g4/g4['NR']['tlive'] + N_er*scale_g4/g4['ER']['tlive'] + N_ng*scale_ng/cap['tlive'])*tlive_PuBe

    ##########
    #Calculate the log probability = log prior + log likelihood
    ll=None
    
    if likelihood=='Norm':
        ll = ll_norm(N_meas[slice(*spec_bounds)],N_pred[slice(*spec_bounds)])
    elif likelihood=='Pois':
        ll = ll_pois(N_meas[slice(*spec_bounds)],N_pred[slice(*spec_bounds)])
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

labels=mcmc_data['labels']
theta_bounds=mcmc_data['theta_bounds']

#Set fit range
E_lim_min=mcmc_data['Efit_min'] #eVee
E_lim_max=mcmc_data['Efit_max'] #eVee
spec_bounds=(np.digitize(E_lim_min,Ebins)-1,np.digitize(E_lim_max,Ebins)-1)
mcmc_data['spec_bounds']=spec_bounds

def Fit_helper(theta):                 
    return calc_log_prob(theta=theta, theta_bounds=theta_bounds,
                         spec_bounds=spec_bounds, likelihood=mcmc_data['likelihood'])

################################################################################
#MCMC Fitting loop
################################################################################

import time
from schwimmbad import MPIPool

nwalkers=mcmc_data['nwalkers']
ndim=mcmc_data['ndim']
nstep=mcmc_data['nstep']

#guesses_s = np.array([0.18, 2e-3, 3.0, 1.0, 1.0]) + np.array([1e-2, 1e-4, 1, 0.1, 0.1]) * np.random.randn(nwalkers, ndim)

#Sample priors uniformly
if mcmc_data['guesses']=='Uniform':
    bounds_low=np.array(theta_bounds)[:,0]
    bounds_hi=np.array(theta_bounds)[:,1]
    mcmc_data['guesses']=(bounds_hi-bounds_low)*np.random.random_sample((nwalkers, ndim))+bounds_low

with MPIPool() as pool:
    if not pool.is_master():
        pool.wait()
        sys.exit(0)
    
    if mcmc_data['moves']=='Default':
        sampler = emcee.EnsembleSampler(nwalkers, ndim, Fit_helper, pool=pool)
    elif mcmc_data['moves']=='DE8020':
        sampler = emcee.EnsembleSampler(nwalkers, ndim, Fit_helper, pool=pool,
                                        moves=[(emcee.moves.DEMove(), 0.8), (emcee.moves.DESnookerMove(), 0.2)])
    else:
        print('Error: mcmc move type not specified')
        
    start = time.time()
    sampler.run_mcmc(mcmc_data['guesses'], nstep, progress=True);
    end = time.time()
    mcmc_data['t_run']=end - start
    
    print('Total run time was {0:.2f} hours.'.format((end - start)/3600))
    
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
    saveFile = open(fname, 'wb')
    
    mcmc_data['sampler']=sampler
    
    pkl.dump(mcmc_data,saveFile)
    saveFile.close()