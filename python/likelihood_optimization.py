def mcmc_data(model=None,update=None) #Data Settings
    mcmc_data = {'g4_load_frac':0.1,
               'cap_load_frac':0.1,
               'cap_sim_file':'../data/v3_400k.pkl',
               'cap_rcapture':0.218,
               ########################## Spectrum Settings ##########################
               'Emax': 2000, #[eVee]
               'Ebins': np.linspace(0,2500,251), #np.linspace(0,2000,201),
               'Efit_min':50, #[eVee]
               'Efit_max':2000, #1750, #[eVee]
               'spectrum_units':'reco-rate', #One of {'counts', 'reco-rate'}
               ########################## Sim Spectra Settings ##########################
               'ER_spec_model':'sim', #One of {'sim', 'flat') to selct from G4 simulation or flat
               'ER_par_labels':[r'$scale_{ER}$'],
               'ER_par_bounds':[(0,20)], #Unitless scaling factor
               'NR_spec_model':'sim', #One of {'sim', 'flat', 'exp') to selct from G4 simulation, flat, or exponential
               'NR_par_labels':[r'$scale_{NR}$'],
               'NR_par_bounds':[(0,20)], #Unitless scaling factor
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
    ########################## Yield Model Settings ##########################
    if model in ['Lind','L','l','lind','lindhard','Lindhard','Lin','lin']:
        mcmc_data.append({
            'Y_model':'Lind',
            'Y_labels': [r'k', r'$F_{NR}$'],
            'Y_bounds': [(0.05,0.3),(0,30)]
        }
    elif model in ['Sor','sor','S','s','Sorensen','sorensen','Sorenson','sorenson']:
        mcmc_data.append({
            'Y_model':'Sor',
            'Y_labels': [r'k', r'q', r'$F_{NR}$'],
            'Y_bounds': [(0.05,0.3),(0,3e-2),(0,30)]
        }
    elif model in ['Chav','chav','Chavaria','chavaria','C','c','Chavarria','chavarria']:
        mcmc_data.append({
            'Y_model':'Chav',
            'Y_labels': [r'k', r'$a^{-1}$', r'$F_{NR}$'],
            'Y_bounds': [(0.05,0.3),(0,2e3),(0,30)]
        }
    elif model in ['AC','ac','A','a']:
        mcmc_data.append({
            'Y_model':'AC',
            'Y_labels': [r'k', r'$\xi$', r'$F_{NR}$'],
            'Y_bounds': [(0.05,0.3),(0,2e3),(0,30)]
        }
    elif !isinstance(model, str):
        raise TypeError('Expected string for model name.',model, type(model))
    
    #Make any other miscellaneous overrides
    if update is not None:
        if !isinstance(update, dict):
            raise TypeError('Expected dictionary for `update` variable.',update,type(update))
        else:
            mcmc_data.update(update)
    
    ################################################################################
    #Parse Options
    ################################################################################
    mcmc_data['labels']=mcmc_data['Y_labels']+mcmc_data['ER_par_labels']+mcmc_data['NR_par_labels']+mcmc_data['NG_par_labels']
    mcmc_data['bounds']=mcmc_data['Y_bounds']+mcmc_data['ER_par_bounds']+mcmc_data['NR_par_bounds']+mcmc_data['NG_par_bounds']
    mcmc_data['ndim']=len(mcmc_data['labels'])
    #Set fit range
    mcmc_data['spec_bounds']=(np.digitize(mcmc_data['Efit_min'],mcmc_data['Ebins'])-1, np.digitize(mcmc_data['Efit_max'],mcmc_data['Ebins'])-1)
            
    ################################################################################
    #Okay we got everything now
    ################################################################################
    return mcmc_data
            
##############################
## Possibly can be embedded ##
##############################
#This stuff is marked as "user doesn't need access; clean up memory use by embedding these in their respective functions.

#Set eVee energy binning
Ebins_ctr=(mcmc_data['Ebins'][:-1]+mcmc_data['Ebins'][1:])/2
#Measured spectra background subtraction
meas=r68.load_measured(verbose=True)
