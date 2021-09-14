#############################
######### Imports ###########
#############################
import R68_spec_tools as Spec
import likelihoods as L
import R68_yield as Yield
import warnings as Warn

#############################
########### Setup ###########
#############################

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
    
    #Make any other miscellaneous overrides before we calculate based on variables.
    if update is not None:
        if !isinstance(update, dict):
            raise TypeError('Expected dictionary for `update` variable.',update,type(update))
        else:
            mcmc_data.update(update)
    
    ################################################################################
    #Parse Options
    ################################################################################
    mcmc_data['labels']=mcmc_data['Y_labels']
            +mcmc_data['ER_par_labels']+mcmc_data['NR_par_labels']+mcmc_data['NG_par_labels']
    mcmc_data['bounds']=mcmc_data['Y_bounds']
            +mcmc_data['ER_par_bounds']+mcmc_data['NR_par_bounds']+mcmc_data['NG_par_bounds']
    mcmc_data['ndim']=len(mcmc_data['labels'])
    #Set fit range
    mcmc_data['spec_bounds']=(np.digitize(mcmc_data['Efit_min'],
                                          mcmc_data['Ebins'])-1,
                              np.digitize(mcmc_data['Efit_max'],mcmc_data['Ebins'])-1)
    #Set eVee energy binning
    mcmc_data['Ebins_ctr']=(mcmc_data['Ebins'][:-1]+mcmc_data['Ebins'][1:])/2
    mcmc_data['N_meas'],mcmc_data['dN_meas']=Spec.doBkgSub(r68.load_measured(), 
                             mcmc_data['Ebins'], mcmc_data['Efit_min'], mcmc_data['Efit_max'],\
                             doEffsyst=mcmc_data['doEffsyst'], doBurstLeaksyst=mcmc_data['doBurstLeaksyst'],\
                             output=mcmc_data['spectrum_units'])
    if mcmc_data['likelihood']=='SNorm':
    #Precalculate split normal likelihood params if we're going to need it
    mcmc_data['SNpars']=L.getSNparsArray(N_meas[slice(*mcmc_data['spec_bounds'])],
                            dN_meas[0][slice(*mcmc_data['spec_bounds'])],
                            dN_meas[1][slice(*mcmc_data['spec_bounds'])]).T
            
    #Load g4 Simulations
    if (mcmc_data['ER_spec_model']=='sim') or (mcmc_data['NR_spec_model']=='sim'):
        mcmc_data['g4']=r68.load_G4(load_frac=mcmc_data['g4_load_frac'])

    #Trim events that won't figure into the analysis range
    #Trimmed sim data
    if (mcmc_data['ER_spec_model']=='sim'):
        mcmc_data['Eee_er']=np.sum(g4['ER']['E'],axis=1)
        mcmc_data['Evec_er_cut']=(Eee_er>10) & (Eee_er<3e3)
        mcmc_data['Evec_er']=g4['ER']['E'][Evec_er_cut]

    if (mcmc_data['NR_spec_model']=='sim'):
        mcmc_data['Eee_nr']=np.sum(g4['NR']['E'],axis=1)
        mcmc_data['Evec_nr_cut']=(Eee_nr>10) & (Eee_nr<30e3)
        mcmc_data['Evec_nr']=g4['NR']['E'][Evec_nr_cut]
            
    mcmc_data['cap']=r68.load_simcap(file=mcmc_data()['cap_sim_file'], 
                    rcapture=mcmc_data()['cap_rcapture'], 
                    load_frac=mcmc_data()['cap_load_frac'])
    ################################################################################
    #Okay we got everything now
    ################################################################################
    return mcmc_data
            
##############################
## Possibly can be embedded ##
##############################
#This stuff is marked as: user doesn't need access; clean up memory use by embedding these in their respective functions.
  
#############################
######### Modeling ##########
#############################
def Y(model):
    tempY = Yield.Yield('Lind',[0.15])
    newmodel = mcmc_data['Y_model'] #give flexibility by "syncincg" variations of string to expected version
    return Yield.Yield(newmodel,np.zeros(tempY.model_npar[newmodel]))

#Version of calc_log_likelihood that is slow because it calls everything each time
#But it requires less setup because it assumes no setup occured.
def calc_log_likelihood_slow(model,theta=[0.2,1,1,1,1],
                        #theta_bounds=tuple(mcmc_data['bounds']), spec_bounds=mcmc_data['spec_bounds'], likelihood=mcmc_data['likelihood'],
                        verbose=False,f_ngzero=False,mcmc_args=None): 
    #Create shorthands for inside of function
    mcmc_data = mcmc_data(model,mcmc_args)
    theta_bounds=tuple(mcmc_data['bounds'])
    spec_bounds=mcmc_data['spec_bounds']
    likelihood=mcmc_data['likelihood']
    
    if mcmc_data['model'] == 'Sor' && len(theta)==5: #Make theta the correct length if left at default
        theta.append(1)
            
            
    global N_pred #Make N_pred callable elsewhere #I should probably change this for this notebook
    N_meas=mcmc_data['N_meas']
    N_er=mcmc_data['N_er']
    Evec_nr=mcmc_data['Evec_nr']
    cap=mcmc_data['cap']
    Y=Y(model)
    
    if f_ngzero:
        theta[-1] = 0
    
    ############
    #Parse fit params
    nYpar=Y.npars

    Y.set_pars(theta[:nYpar])
    F_NR=theta[nYpar] #theta[2] 
    if verbose:
        print(Y.pars,F_NR)
    
    if not Y.solve():
        return -np.inf
        
    #Grab correct ER scaling factor from theta
    if mcmc_data['ER_spec_model']=='sim':
        if '$scale_{G4}$' in mcmc_data['labels']:
            scale_er=theta[mcmc_data['labels'].index('$scale_{G4}$')]/g4['ER']['tlive']
        else:
            scale_er=theta[mcmc_data['labels'].index('$scale_{ER}$')]/g4['ER']['tlive']
        
    ##########
    #Build the spectra    
    if mcmc_data['spectrum_units']=='reco-rate':
        #Don't apply any efficiency effects to simulated spectrum
        #NR
        if mcmc_data['NR_spec_model']=='sim':
            N_nr=spec.buildAvgSimSpectrum_ee(Ebins=mcmc_data['Ebins'], 
                                             Evec=Evec_nr, Yield=Y, F=F_NR, scale=1,\
                                             doDetRes=mcmc_data['doDetRes'], fpeak=mcmc_data['fpeak'],\
                                             doEffs=False)
            if verbose:
                print("\nN_nr INPUTS")
                print("-------------")
                print("Ebins:",mcmc_data['Ebins'])
                print("Evec_nr:",Evec_nr)
                print("Yield:",Y.pars)
                print("F_NR:",F_NR)
                print("doDetRes:",mcmc_data['doDetRes'])
                print("fpeak:",mcmc_data['fpeak'],'\n')

            if '$scale_{G4}$' in mcmc_data['labels']:
                scale_nr=theta[mcmc_data['labels'].index('$scale_{G4}$')]/g4['NR']['tlive']
            else:
                scale_nr=theta[mcmc_data['labels'].index('$scale_{NR}$')]/g4['NR']['tlive']
        elif mcmc_data['NR_spec_model']=='exp':
            N_nr=np.exp(-Ebins_ctr/theta[mcmc_data['labels'].index(r'$E0_{NR}$')])
            scale_nr=theta[mcmc_data['labels'].index(r'$R0_{NR}$')]
        
        if verbose:
            print('spec_units == reco-rate\n')
            print('N_ng INPUTS')
            print('-----------')
            print('Ebins:',mcmc_data['Ebins'])
            print('Evec:',cap['E'])
            print('dEvec:',cap['dE'])
            print('Yield:',Y.pars)
            print('F:',F_NR)
            print('doDetRes:',mcmc_data['doDetRes'])
            print('fpeak:',mcmc_data['fpeak'],'\n')
        #(n,gamma)
        N_ng=spec.buildAvgSimSpectrum_ee_composite(Ebins=mcmc_data['Ebins'], Evec=cap['E'], dEvec=cap['dE'],\
                                                   Yield=Y, F=F_NR, scale=1, doDetRes=mcmc_data['doDetRes'],\
                                                   fpeak=mcmc_data['fpeak'], doEffs=False,verbose=verbose)
        scale_ng=theta[mcmc_data['labels'].index('$scale_{ng}$')]/cap['tlive']
        
        #Calculate rate (though we'll still call them N_* just to be confusing)
        N_pred = (N_nr*scale_nr + 
                  N_er*scale_er + 
                  N_ng*scale_ng)
        if verbose:
            print("N_ng:",N_ng[0:3]) 
            print("N_nr:", N_nr[0:3])
            print("F_NR:", F_NR)
            print("N_pred:",N_pred[0])

    ##########
    #Calculate the log probability = log prior + log likelihood
    ll=None
    
    if likelihood=='Norm':
        ll = ll_norm(N_meas[slice(*spec_bounds)],N_pred[slice(*spec_bounds)])
    elif likelihood=='Pois':
        ll = ll_pois(N_meas[slice(*spec_bounds)],N_pred[slice(*spec_bounds)])
    elif likelihood=='SNorm':
        #Precalculate split normal likelihood params if we're going to need it
        SNpars=L.getSNparsArray(N_meas[slice(*spec_bounds)],
                                dN_meas[0][slice(*spec_bounds)],
                                dN_meas[1][slice(*spec_bounds)])
        SNpars=SNpars.T
        ll = L.ll_SNorm(N_pred[slice(*spec_bounds)],*SNpars) #4 parameters
        if verbose:
            print("SNpars.T:",SNpars[:,-5:])
            print("N_pred:",N_pred[slice(*spec_bounds)][:5])
    else:
        print('Error: Bad likelihood')
        return None
    
    if not np.isfinite(ll):
        return -np.inf
    
    return ll #When optimizing, I only want to optimize the log likelihood, not of prior*likelihood.
    #So I've slightly altered this function. Also shrunk down to SNorm for ease of troubleshooting

            
            
#"Normal speed" version of calc_log_likelihood that requires mcmc_data to exist
def calc_log_likelihood(model,theta=[0.2,1,1,1,1],
                        #theta_bounds=tuple(mcmc_data['bounds']), spec_bounds=mcmc_data['spec_bounds'], likelihood=mcmc_data['likelihood'],
                        verbose=False,f_ngzero=False): 
    #Create shorthands for inside of function
    #mcmc_data = mcmc_data(model,mcmc_args)
    theta_bounds=tuple(mcmc_data['bounds'])
    spec_bounds=mcmc_data['spec_bounds']
    likelihood=mcmc_data['likelihood']
    
    if mcmc_data['model'] == 'Sor' && len(theta)==5: #Make theta the correct length if left at default
        theta.append(1)
            
            
    global N_pred #Make N_pred callable elsewhere #I should probably change this for this notebook
    N_meas=mcmc_data['N_meas']
    N_er=mcmc_data['N_er']
    Evec_nr=mcmc_data['Evec_nr']
    cap=mcmc_data['cap']
    Y=Y(model)
    
    if f_ngzero:
        theta[-1] = 0
    
    ############
    #Parse fit params
    nYpar=Y.npars

    Y.set_pars(theta[:nYpar])
    F_NR=theta[nYpar] #theta[2] 
    if verbose:
        print(Y.pars,F_NR)
    
    if not Y.solve():
        return -np.inf
        
    #Grab correct ER scaling factor from theta
    if mcmc_data['ER_spec_model']=='sim':
        if '$scale_{G4}$' in mcmc_data['labels']:
            scale_er=theta[mcmc_data['labels'].index('$scale_{G4}$')]/g4['ER']['tlive']
        else:
            scale_er=theta[mcmc_data['labels'].index('$scale_{ER}$')]/g4['ER']['tlive']
        
    ##########
    #Build the spectra    
    if mcmc_data['spectrum_units']=='reco-rate':
        #Don't apply any efficiency effects to simulated spectrum
        #NR
        if mcmc_data['NR_spec_model']=='sim':
            N_nr=spec.buildAvgSimSpectrum_ee(Ebins=mcmc_data['Ebins'], 
                                             Evec=Evec_nr, Yield=Y, F=F_NR, scale=1,\
                                             doDetRes=mcmc_data['doDetRes'], fpeak=mcmc_data['fpeak'],\
                                             doEffs=False)
            if verbose:
                print("\nN_nr INPUTS")
                print("-------------")
                print("Ebins:",mcmc_data['Ebins'])
                print("Evec_nr:",Evec_nr)
                print("Yield:",Y.pars)
                print("F_NR:",F_NR)
                print("doDetRes:",mcmc_data['doDetRes'])
                print("fpeak:",mcmc_data['fpeak'],'\n')

            if '$scale_{G4}$' in mcmc_data['labels']:
                scale_nr=theta[mcmc_data['labels'].index('$scale_{G4}$')]/g4['NR']['tlive']
            else:
                scale_nr=theta[mcmc_data['labels'].index('$scale_{NR}$')]/g4['NR']['tlive']
        elif mcmc_data['NR_spec_model']=='exp':
            N_nr=np.exp(-Ebins_ctr/theta[mcmc_data['labels'].index(r'$E0_{NR}$')])
            scale_nr=theta[mcmc_data['labels'].index(r'$R0_{NR}$')]
        
        if verbose:
            print('spec_units == reco-rate\n')
            print('N_ng INPUTS')
            print('-----------')
            print('Ebins:',mcmc_data['Ebins'])
            print('Evec:',cap['E'])
            print('dEvec:',cap['dE'])
            print('Yield:',Y.pars)
            print('F:',F_NR)
            print('doDetRes:',mcmc_data['doDetRes'])
            print('fpeak:',mcmc_data['fpeak'],'\n')
        #(n,gamma)
        N_ng=spec.buildAvgSimSpectrum_ee_composite(Ebins=mcmc_data['Ebins'], Evec=cap['E'], dEvec=cap['dE'],\
                                                   Yield=Y, F=F_NR, scale=1, doDetRes=mcmc_data['doDetRes'],\
                                                   fpeak=mcmc_data['fpeak'], doEffs=False,verbose=verbose)
        scale_ng=theta[mcmc_data['labels'].index('$scale_{ng}$')]/cap['tlive']
        
        #Calculate rate (though we'll still call them N_* just to be confusing)
        N_pred = (N_nr*scale_nr + 
                  N_er*scale_er + 
                  N_ng*scale_ng)
        if verbose:
            print("N_ng:",N_ng[0:3]) 
            print("N_nr:", N_nr[0:3])
            print("F_NR:", F_NR)
            print("N_pred:",N_pred[0])

    ##########
    #Calculate the log probability = log prior + log likelihood
    ll=None
    
    if likelihood=='Norm':
        ll = ll_norm(N_meas[slice(*spec_bounds)],N_pred[slice(*spec_bounds)])
    elif likelihood=='Pois':
        ll = ll_pois(N_meas[slice(*spec_bounds)],N_pred[slice(*spec_bounds)])
    elif likelihood=='SNorm':
        #Precalculate split normal likelihood params if we're going to need it
        SNpars=L.getSNparsArray(N_meas[slice(*spec_bounds)],
                                dN_meas[0][slice(*spec_bounds)],
                                dN_meas[1][slice(*spec_bounds)])
        SNpars=SNpars.T
        ll = L.ll_SNorm(N_pred[slice(*spec_bounds)],*SNpars) #4 parameters
        if verbose:
            print("SNpars.T:",SNpars[:,-5:])
            print("N_pred:",N_pred[slice(*spec_bounds)][:5])
    else:
        print('Error: Bad likelihood')
        return None
    
    if not np.isfinite(ll):
        return -np.inf
    
    return ll #When optimizing, I only want to optimize the log likelihood, not of prior*likelihood.
    #So I've slightly altered this function. Also shrunk down to SNorm for ease of troubleshooting
            
def optll(model,initial,bounds,mcmc_args=None):
    mcmc_data = mcmc_data(model,mcmc_args)
    nll = lambda *args: -calc_log_likelihood(model,mcmc_args=mcmc_args,*args)
    soln = minimize(nll,initial,method='SLSQP',bounds=bounds)
    params = {} #Initialize empty dictionary
    if mcmc_data['model'] == 'Sor':
        for i,x in enumerate('k','q','F_NR','f_ER','f_NR','f_ng'):
            params.update(x:params[i])
    elif mcmc_data['model'] == 'Lind':
        for i,x in enumerate('k','F_NR','f_ER','f_NR','f_ng'):
            params.update(x:params[i])
    elif mcmc_data['model'] == 'Chav':
        Warn.warn("Individual variables for model 'Chav' not yet implemented. Retrieve from list (`optll()[0][i]`).")
    elif mcmc_data['model'] == 'AC':
        Warn.warn("Individual variables for model 'AC' not yet implemented. Retrieve from list (`optll()[0][i]`).")
    return soln,params

def print_stats(model,params,comment=None,mcmc_args=None):
    mcmc_data=mcmc_data(model,mcmc_args)
    model = mcmc_data['model']
    print("MODEL:",model,"\n")

    print("STATISTICS")
    print("------------")
    print("Log likelihood:",calc_log_likelihood(model,theta=params))

    if mcmc_data['likelihood'] == 'SNorm':
        chi = L.chisq(N_meas[slice(*mcmc_data['spec_bounds'])],N_pred[slice(*mcmc_data['spec_bounds'])],0.5*(dN_meas[0]+dN_meas[1])[slice(*mcmc_data['spec_bounds'])])
    else:
        raise ValueError('Available likelihood models: SNorm.',mcmc_data['likelihood'])
    dof=np.diff(mcmc_data['spec_bounds'])[0]-mcmc_data['ndim']

    print('Chisq:',chi)
    print('Chisq/DoF:',chi/dof,"\n")

    print("PARAMETERS")
    print("------------")
    if model == 'Sor':
        k,q,F_NR,f_ER,f_NR,f_ng = *params
        print("k = {0:.3f}".format(k))
        print("q = {0:.5e}".format(q))
        print("F_NR = {0:.3f}".format(F_NR))
        print("f_ER = {0:.3f}".format(f_ER))
        print("f_NR = {0:.3f}".format(f_NR))
        print("f_ng = {0:.3f}\n".format(f_ng))
    elif model == 'Lind':
        k,F_NR,f_ER,f_NR,f_ng = *params
        print("k = {0:.3f}".format(k))
        print("F_NR = {0:.3f}".format(F_NR))
        print("f_ER = {0:.3f}".format(f_ER))
        print("f_NR = {0:.3f}".format(f_NR))
        print("f_ng = {0:.3f}\n".format(f_ng))
    elif model == 'Chav':
        Warn.warn("Parameter printing for model 'Chav' not yet implemented.")
    elif model == 'AC':
        Warn.warn("Parameter printing for model 'AC' not yet implemented.")
    
    if comment != None:
        if !isinstance(comment, str):
            raise TypeError('Expected string for comment.',model, type(model))
            
def get_stats(model,initial,bounds,mcmc_args=None,comment=None):
    mcmc_data = mcmc_data(model,mcmc_args)
    params = optll(model,initial,bounds,mcmc_args)[0]
    print_stats(model,params,comment,mcmc_args)