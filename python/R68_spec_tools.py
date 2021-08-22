###########################################################################
#Methods for working with spectra in R68
###########################################################################
import constants as const
import R68_efficiencies as eff
from scipy.special import erf
import numpy as np
import R68_load

###########################################################################
#Do background subtraction of measured data
#
#meas: dictionary of measured events. Should be the output of R68_load.load_measured
#Ebins: [eVee] Histogram energy bins. 1-d, length M+1
#Efit_min: [eVee] Minimum energy over which the fit will be performed (min analysis threshold)
#Efit_max: [eVee] Maximum energy over which the fit will be performed (max analysis threshold)
#doEffsyst: Include efficiency uncertainties (from write and cut efficiencies)
#doBurstLeaksyst: Include systematic uncertainty of leakage from burst cut
#output: What quantities to output. One of {'counts', 'reco-rate'}. 
#  'counts' will return the background-subtracted PuBe counts: N_PuBe-N_Bkg*ratio_Tlive*ratio_write*ratio_cutEff
#  'reco-rate' will return the reconstructed measured rate: 
#      N_PuBe/(Tlive_PuBe*eff_{write,PuBe}*eff_{cut,PuBe}*eff_trig)  -
#        N_Bkg/(Tlive_Bkg*eff_{write,Bkg}*eff_{cut,Bkg}*eff_trig)
#
#Returns: (N_meas, dN_meas) The background-subtracted measured spectrum and its uncertainty.
# Units are either counts/bin or rate.
# N_meas is bin counts. 1-d, length M
# dN_meas is (high,low) uncertainties, 2-d, each length M.
# Uncertatinty always includes Poisson uncertainties of Bkg and PuBe count rates at least. 
# More effects can be enabled with the doEffsyst and doLERsyst tags
###########################################################################
def doBkgSub(meas, Ebins, Efit_min, Efit_max, doEffsyst=True, doBurstLeaksyst=False, output='counts'):
    Ebins_ctr=(Ebins[:-1]+Ebins[1:])/2
    
    N_meas_PuBe,_ = np.histogram(meas['PuBe']['E'],bins=Ebins)
    N_meas_Bkg,_ = np.histogram(meas['Bkg']['E'],bins=Ebins)
    #Count uncertainties are Poisson
    dN_meas_PuBe_Pois=np.sqrt(N_meas_PuBe)
    dN_meas_Bkg_Pois=np.sqrt(N_meas_Bkg)

    #Call the factor of livetime*efficiencies TE_PuBe and TE_bkg
    #We either want to return 
    #R_meas = N_PuBe/TE_PuBe - N_Bkg/TE_Bkg
    # or
    #N_meas = N_PuBe-N_Bkg*TE_PuBe/TE_Bkg
    
    TE_PuBe = meas['PuBe']['tlive']*eff.eff_write*eff.cutEffFit(Ebins_ctr)*eff.trigEff(Ebins_ctr)
    TE_Bkg = meas['Bkg']['tlive']*eff.eff_write_bkg*eff.cutEffFit_bkg(Ebins_ctr)*eff.trigEff(Ebins_ctr)
    
    dTE_PuBe = TE_PuBe*np.sqrt( (eff.deff_write/eff.eff_write)**2 +\
                               (eff.dcutEffFit(Ebins_ctr)/eff.cutEffFit(Ebins_ctr))**2 +\
                               (eff.dtrigEff(Ebins_ctr)/eff.trigEff(Ebins_ctr))**2 )
    
    dTE_Bkg = TE_PuBe*np.sqrt( (eff.deff_write_bkg/eff.eff_write_bkg)**2 +\
                               (eff.dcutEffFit_bkg(Ebins_ctr)/eff.cutEffFit_bkg(Ebins_ctr))**2 +\
                               (eff.dtrigEff(Ebins_ctr)/eff.trigEff(Ebins_ctr))**2 )
    if output=='reco-rate':
        #Reconstructed rate
        R_meas_PuBe = N_meas_PuBe/TE_PuBe
        R_meas_Bkg = N_meas_Bkg/TE_Bkg
        
        dR_meas_PuBe_hi = R_meas_PuBe*np.sqrt( (dN_meas_PuBe_Pois/N_meas_PuBe)**2 +\
                                              doEffsyst*(dTE_PuBe/TE_PuBe)**2 )
        dR_meas_PuBe_low = R_meas_PuBe*np.sqrt( (dN_meas_PuBe_Pois/N_meas_PuBe)**2 +\
                                               doBurstLeaksyst*(eff.dtrigburstLeak(Ebins_ctr)/N_meas_PuBe)**2 +\
                                               doEffsyst*(dTE_PuBe/TE_PuBe)**2 )
        
        dR_meas_Bkg = R_meas_Bkg*np.sqrt( (dN_meas_Bkg_Pois/N_meas_Bkg)**2 +\
                                         doEffsyst*(dTE_Bkg/TE_Bkg)**2 )

        R_meas = R_meas_PuBe - R_meas_Bkg
        
        dR_meas_hi = np.sqrt(dR_meas_PuBe_hi**2 + dR_meas_Bkg**2)
        dR_meas_low = np.sqrt(dR_meas_PuBe_low**2 + dR_meas_Bkg**2)
        
        return (R_meas, np.array([dR_meas_hi, dR_meas_low]))
    
    elif output=='counts':
        #Bkg-subtracted measured PuBe signal
        N_bkg_scaled=N_meas_Bkg*TE_PuBe/TE_Bkg
        
        #Here is a little kludge
        #We don't have a way to include the efficiency uncertainties for the simulated specrtum,
        #  so we'll fold those effects into the measured spectra here. This isn't quite right, but
        #  should help.
        #Include uncertainty from efficiencies
        dN_meas_PuBe_hi = N_meas_PuBe*np.sqrt( (dN_meas_PuBe_Pois/N_meas_PuBe)**2 +\
                                           doEffsyst*(eff.deff_write/eff.eff_write)**2 +\
                                           doEffsyst*(eff.dcutEffFit(Ebins_ctr)/eff.cutEffFit(Ebins_ctr))**2)
        dN_meas_PuBe_low = N_meas_PuBe*np.sqrt( (dN_meas_PuBe_Pois/N_meas_PuBe)**2 +\
                                           doBurstLeaksyst*(eff.dtrigburstLeak(Ebins_ctr)/N_meas_PuBe)**2 +\
                                           doEffsyst*(eff.deff_write/eff.eff_write)**2 +\
                                           doEffsyst*(eff.dcutEffFit(Ebins_ctr)/eff.cutEffFit(Ebins_ctr))**2)
        dN_meas_Bkg = N_meas_Bkg*np.sqrt( (dN_meas_Bkg_Pois/N_meas_Bkg)**2 +\
                                         doEffsyst*(eff.deff_write_bkg/eff.eff_write_bkg)**2 +\
                                         doEffsyst*(eff.dcutEffFit_bkg(Ebins_ctr)/eff.cutEffFit_bkg(Ebins_ctr))**2)
        
        dN_bkg_scaled=N_bkg_scaled*np.sqrt( (dN_meas_Bkg/N_meas_Bkg)**2 +\
                                          doEffsyst*(dTE_PuBe/TE_PuBe)**2 +\
                                          doEffsyst*(dTE_Bkg/TE_Bkg)**2 )
        
        N_meas = N_meas_PuBe - N_bkg_scaled
        dN_meas_hi = np.sqrt( dN_meas_PuBe_hi**2 + dN_bkg_scaled**2 )
        dN_meas_low = np.sqrt( dN_meas_PuBe_low**2 + dN_bkg_scaled**2 )

        return (N_meas, np.array([dN_meas_hi, dN_meas_low]))
        
    else:
        print(f'Error in R68_spec_tools.doBkgSub: bad "output" keyword: {output}')
        print('Should be one of: {"counts", "reco-rate"}')
        return None
###########################################################################
#Apply yield, resolution, efficiencies etc. to simulated data to get an instance of a simulated E_ee spectrum

#Ebins: Bin edges in [eV]
#Evec_nr: Array of NR hit vectors
#Evec_er: Array of ER hit vectors
#Evec_ng: Array of (n,gamma) step energies
#dEvec_ng: Array of (n,gamma) step deposited energies
#Yield: Yield object with a .calc method as defined in R68_yield.py
#F_NR: Nuclear recoil Fano factor
#scale_g4: Scaling of NR+ER sim spectrum 
#scale_ng(n,gamma) spectrum scaling factor
#doDetRes: whether to apply detector resolution function
#seed: random number generator seed
#
#Returns: (n_Eee_nr, n_Eee_er, n_Eee_ng) The binned spectra of NR, ER, and (n,gamma) events in units of [counts/bin]

############################################################################################################
############# DEFUNCT! Use buildAvgSimSpectrum_ee or buildAvgSimSpectrum_ee_composite instead! #############
############################################################################################################
def buildSimSpectra_ee(Ebins, Evec_nr, Evec_er, Evec_ng, dEvec_ng, Yield, F_NR, scale_g4=1, scale_ng=1, doDetRes=True, seed=1):
    
    V=const.V
    G_NTL=const.G_NTL
    eps=const.eps
    F=const.F
    
    np.random.seed(seed)

    ###############
    #Yield, Fano, Resolution
    ###############
    
    ###############
    #NR, No capture
    Evec_nr_eion = Evec_nr*Yield.calc(Evec_nr) #Ionization energy from NRs
    #Number of e/h pairs produced in each event
    Nvec_nr = getNeh(Evec_nr_eion,eps,F_NR)
    #get the total phonon energies
    Etvec_nr = Evec_nr + Nvec_nr*V
    #Sum event and convert to eVee scale by making the gamma assumption
    Eee_nr  = np.sum(Etvec_nr,1)/G_NTL
    #Apply trigger logic to randomly drop events according to their energies
    Eee_nr = Eee_nr[passTrig(Eee_nr)]
    #Detector resolution
    if doDetRes:
        Eee_nr = getSmeared(Eee_nr)
        
    ###############
    #ER, No capture
    E_er_eion = np.sum(Evec_er,1) #Ionization energy from ERs
    #Number of e/h pairs produced in event
    N_er = getNeh(E_er_eion,eps,F)
    #get the total phonon energy
    Et_er = N_er*V + E_er_eion
    #convert to eVee scale
    Eee_er  = Et_er/G_NTL
    #Apply trigger logic to randomly drop events according to their energies
    Eee_er = Eee_er[passTrig(Eee_er)]
    #Detector resolution
    if doDetRes:
        Eee_er = getSmeared(Eee_er)
    
    ###############
    #(n,gamma)
    #Ionization energy available in each step
    Evec_ng_eion = Evec_ng*Yield.calc(Evec_ng) - (Evec_ng-dEvec_ng)*Yield.calc(Evec_ng-dEvec_ng)
    #Number of e/h pairs produced in each step
    Nvec_ng = getNeh(Evec_ng_eion,eps,F_NR)
    #Total phonon energy
    Etvec_ng = Nvec_ng*V + dEvec_ng
    Et_ng = np.sum(Etvec_ng,1)
    #convert to eVee scale by making the gamma assumption
    Eee_ng  = Et_ng/G_NTL
    #Apply trigger logic to randomly drop events according to their energies
    Eee_ng = Eee_ng[passTrig(Eee_ng)]
    #Detector resolution
    if doDetRes:
        Eee_ng = getSmeared(Eee_ng)

    ###############
    #Binning
    ###############
    Ebin_ctr=(Ebins[:-1]+Ebins[1:])/2
    
    n_Eee_nr,_ = np.histogram(Eee_nr,bins=Ebins) #[counts/bin]
    n_Eee_er,_ = np.histogram(Eee_er,bins=Ebins)
    n_Eee_ng,_ = np.histogram(Eee_ng,bins=Ebins)

    ###############
    #Efficiencies
    ###############
    
    #Analysis cut efficiencies
    eff_cuts = eff.cutEff(Ebin_ctr)
    
    n_Eee_nr = n_Eee_nr*eff_cuts*scale_g4
    n_Eee_er = n_Eee_er*eff_cuts*scale_g4
    n_Eee_ng = n_Eee_ng*eff_cuts*scale_ng
    
    return (n_Eee_nr, n_Eee_er, n_Eee_ng)


###########################################################################
#Apply yield, resolution, efficiencies etc. to simulated data to get an average spectrum
#
# This version convolves fano and resolution functions with the spectrum, as opposed to drawing invividual values from 
# probability distributions for each hit.
# As such, this function returns the average spectra, which will not, in general, have integer counts
#
#Ebins: Bin edges in [eV]
#Evec_nr: Array of NR hit vectors
#Evec_er: Array of ER hit vectors
#Evec_ng: Array of (n,gamma) step energies
#dEvec_ng: Array of (n,gamma) step deposited energies
#Yield: Yield object with a .calc method as defined in R68_yield.py
#F_NR: Nuclear recoil Fano factor
#scale_g4: Scaling of NR+ER sim spectrum 
#scale_ng(n,gamma) spectrum scaling factor
#doDetRes: whether to apply detector resolution function
#
#Returns: (n_Eee_nr, n_Eee_er, n_Eee_ng) The binned spectra of NR, ER, and (n,gamma) events in units of [counts/bin]
#
#TODO: -Handle low energy Neh consistent with model used in buildSpectra
############################################################################################################
############# DEFUNCT! Use buildAvgSimSpectrum_ee or buildAvgSimSpectrum_ee_composite instead! #############
############################################################################################################
def buildAvgSimSpectra_ee(Ebins, Evec_nr, Evec_er, Evec_ng, dEvec_ng, Yield, F_NR, scale_g4=1, scale_ng=1, doDetRes=True, fpeak=0.753):

    V=const.V
    G_NTL=const.G_NTL
    eps=const.eps
    F=const.F
    
    #Params from Matt's Bkg resolution fit:
    #https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_panda:calibration#resolution_versus_energy
    sigma0=const.sigma0 #eV 
    B=const.B #This includes FANO
    #B_1=B-F*eps #Actually only true for V->inf (http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/Peak_Widths.pdf)
    B_1=B-F*V**2/eps/(1+V/eps)**2
    A=const.A

    ###############
    #Binning
    ###############
    Ebin_ctr = (Ebins[:-1]+Ebins[1:])/2

    ###############
    #Yield, Fano
    ###############

    ################
    #NR, No capture
    #<E_ee> for each hit
    EeeVec_nr = Evec_nr*(1 + Yield.calc(Evec_nr)*V/eps) / G_NTL
    
    #Fano width in Evee for each hit
    SigmaEeeSqVec_nr = Evec_nr*Yield.calc(Evec_nr)*(V**2)*F_NR/eps/(G_NTL**2)
    
    #energy==0 (no hit) may still generate Yield!=0, take care of that here
    SigmaEeeSqVec_nr = np.where(Evec_nr>0, SigmaEeeSqVec_nr, 0)
    
    #Sum hits assuming they're uncorrelated
    Eee_nr=np.sum(EeeVec_nr,1)
    SigmaEee_nr = np.sqrt(np.sum(SigmaEeeSqVec_nr,1)) 
    
    #Reshape so we can broadcast against bin centers
    Eee_nr = Eee_nr[:,np.newaxis]
    SigmaEee_nr =SigmaEee_nr[:,np.newaxis]
    
    #Weighted spectra of Eee energies measured
    #nvec_Eee_nr = gausInt(Eee_nr, SigmaEee_nr, Ebins[:-1], Ebins[1:])
    nvec_Eee_nr = gausInt_fast(Eee_nr, SigmaEee_nr, Ebins)
    #NR spectrum
    n_Eee_nr = np.sum(nvec_Eee_nr,0)

    
    ################
    #ER, No capture

    #Recoil energy, Er, in each event
    Eee_er = np.sum(Evec_er,1)
    #Fano width in eVee
    SigmaEee_er = np.sqrt(Eee_er*(V**2)*F/eps/(G_NTL**2))

    #Reshape so we can broadcast against bin centers
    Eee_er = Eee_er[:,np.newaxis]
    SigmaEee_er = SigmaEee_er[:,np.newaxis]
    
    #Weighted spectra of Eee energies
    #nvec_Eee_er = gausInt(Eee_er, SigmaEee_er, Ebins[:-1], Ebins[1:])
    nvec_Eee_er = gausInt_fast(Eee_er, SigmaEee_er, Ebins)
    #ER spectrum
    n_Eee_er = np.sum(nvec_Eee_er,0)
    
    ################
    #N,gamma
    #Convert Simulated n,gamma hits to eVee events
    
    #Ionization energy available in each step
    Evec_ng_eion = Evec_ng*Yield.calc(Evec_ng) - (Evec_ng-dEvec_ng)*Yield.calc(Evec_ng-dEvec_ng)
    #Electon-equivalent
    EeeVec_ng = (dEvec_ng + Evec_ng_eion*V/eps) / G_NTL
    #Fano resolution in eVee
    #TODO: This is not quite correct if dE/E<<1, but may be close enough since that's rare?
    SigmaEeeSqVec_ng = Evec_ng_eion*(V**2)*F_NR/eps/(G_NTL**2)

    #energy==0 (no hit) may still generate Yield!=0, take care of that here
    SigmaEeeSqVec_ng = np.where(Evec_ng_eion>0, SigmaEeeSqVec_ng, 0)

    #Sum hits assuming they're uncorrelated
    Eee_ng = np.sum(EeeVec_ng,1)
    SigmaEee_ng = np.sqrt(np.sum(SigmaEeeSqVec_ng,1))  
    
    #Reshape so we can broadcast against bin centers
    Eee_ng = Eee_ng[:,np.newaxis]
    SigmaEee_ng = SigmaEee_ng[:,np.newaxis]

    #Weighted spectra of Eee energies measured
    #nvec_Eee_ng = gausInt(Eee_ng, SigmaEee_ng, Ebins[:-1], Ebins[1:])
    nvec_Eee_ng = gausInt_fast(Eee_ng, SigmaEee_ng, Ebins)
    #(n,gamma) spectrum
    n_Eee_ng = np.sum(nvec_Eee_ng,0)

    ###############
    #Vbias Smearing
    ###############
    #Currently only for studying systematic uncertainties, smearing model is an over-estimation
    #See EPot_Model.ipynb
    if fpeak<1.0:
        expdelta_kernel=expdelta_int(Ebin_ctr[:,np.newaxis],fpeak,4.20823831,Ebins)
    
        #Smear spectrum using resolution function
        #NR 
        n_Eee_nr_smeared_vec = n_Eee_nr[:,np.newaxis]*expdelta_kernel
        n_Eee_nr = np.sum(n_Eee_nr_smeared_vec,0)

        #ER
        n_Eee_er_smeared_vec = n_Eee_er[:,np.newaxis]*expdelta_kernel
        n_Eee_er = np.sum(n_Eee_er_smeared_vec,0)

        #(n,gamma)
        n_Eee_ng_smeared_vec = n_Eee_ng[:,np.newaxis]*expdelta_kernel
        n_Eee_ng = np.sum(n_Eee_ng_smeared_vec,0)
    
    ###############
    #Trigger Efficiency
    ###############
    tEff=eff.trigEff(Ebin_ctr)
    n_Eee_nr = n_Eee_nr*tEff
    n_Eee_er = n_Eee_er*tEff
    n_Eee_ng  = n_Eee_ng*tEff
    
    ###############
    #Detector Resolution
    ###############
    if doDetRes:
        #Calc gaus kernels once and reuse them for each broadcast
        #gaus_kernel = gausInt(Ebin_ctr[:,np.newaxis], sigma_ee(Ebin_ctr[:,np.newaxis],sigma0,B_1,A), Ebins[:-1], Ebins[1:])
        gaus_kernel = gausInt_fast(Ebin_ctr[:,np.newaxis], sigma_ee(Ebin_ctr[:,np.newaxis],sigma0,B_1,A), Ebins)
        
        #Smear spectrum using resolution function
        #NR 
        n_Eee_nr_smeared_vec = n_Eee_nr[:,np.newaxis]*gaus_kernel
        n_Eee_nr = np.sum(n_Eee_nr_smeared_vec,0)
        
        #ER
        n_Eee_er_smeared_vec = n_Eee_er[:,np.newaxis]*gaus_kernel
        n_Eee_er = np.sum(n_Eee_er_smeared_vec,0)

        #(n,gamma)
        n_Eee_ng_smeared_vec = n_Eee_ng[:,np.newaxis]*gaus_kernel
        n_Eee_ng = np.sum(n_Eee_ng_smeared_vec,0)
    
    ###############
    #Efficiencies
    ###############
    
    #Analysis cut efficiencies
    eff_cuts = eff.cutEff(Ebin_ctr)
    
    n_Eee_nr = n_Eee_nr*eff_cuts*scale_g4
    n_Eee_er = n_Eee_er*eff_cuts*scale_g4
    n_Eee_ng = n_Eee_ng*eff_cuts*scale_ng
    
    return (n_Eee_nr, n_Eee_er, n_Eee_ng)


###########################################################################
#Apply yield, resolution, efficiencies etc. to simulated ER or NR hit data to get an average spectrum
##########################
#NOTE: For composite NR's like (n,gamma), use buildAvgSimSpectrum_ee_composite instead!
##########################
#
#This version just takes a single data type (either ER or NR) at a time.
#
#Ebins: Bin edges in [eV]
#Evec: Array of hit vectors
#Yield: Either a constant or a Yield object with a .calc method as defined in R68_yield.py
#F: Fano factor
#scale: Scaling of spectrum 
#doDetRes: whether to apply detector resolution function
#fpeak: Vbias smearing effect parameter. Fraction of events which are NOT smeared. (1.0=no smearing) See EPot_Model.ipynb
#doEffs: Apply trigger and cut efficiency curves
#Returns: n_Eee: The binned spectra of events in units of [counts/bin]
#
def buildAvgSimSpectrum_ee(Ebins, Evec, Yield, F, scale, doDetRes=True, fpeak=1.0, doEffs=True):

    V=const.V
    G_NTL=const.G_NTL
    eps=const.eps
    F_ER=const.F
    
    #Params from Matt's Bkg resolution fit:
    #https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_panda:calibration#resolution_versus_energy
    sigma0=const.sigma0 #eV 
    B=const.B #This includes FANO
    #B_1=B-F*eps #Actually only true for V->inf (http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/Peak_Widths.pdf)
    B_1=B-F_ER*V**2/eps/(1+V/eps)**2
    A=const.A

    ###############
    #Binning
    ###############
    Ebin_ctr = (Ebins[:-1]+Ebins[1:])/2

    ###############
    #Yield, Fano
    ###############

    if isinstance(Yield,(float,int)):
        #<E_ee> for each hit
        EeeVec = Evec*(1 + Yield*V/eps) / G_NTL
        #Fano width in Evee^2 for each hit
        SigmaEeeSqVec = Evec*Yield*(V**2)*F/eps/(G_NTL**2)
    else:
        #<E_ee> for each hit
        EeeVec = Evec*(1 + Yield.calc(Evec)*V/eps) / G_NTL
        #Fano width in Evee^2 for each hit
        SigmaEeeSqVec = Evec*Yield.calc(Evec)*(V**2)*F/eps/(G_NTL**2)
    
    
    #energy==0 (no hit) may still generate Yield!=0, take care of that here
    SigmaEeeSqVec[Evec<=0]=0.0
    
    #Sum hits assuming they're uncorrelated and reshape so we can broadcast against bin centers
    Eee=np.sum(EeeVec,1)[:,np.newaxis]
    SigmaEee = np.sqrt(np.sum(SigmaEeeSqVec,1))[:,np.newaxis]
    
    #Weighted spectra of Eee energies measured
    nvec_Eee = gausInt_fast(Eee, SigmaEee, Ebins)
    #Spectrum
    n_Eee = np.sum(nvec_Eee,0)

        
    ###############
    #Vbias Smearing
    ###############
    #Currently only for studying systematic uncertainties, smearing model is an over-estimation
    #See EPot_Model.ipynb
    if fpeak<1.0:
        #Smear spectrum using resolution function
        expdelta_kernel=expdelta_int(Ebin_ctr[:,np.newaxis],fpeak,4.20823831,Ebins)
        n_Eee = np.sum(n_Eee[:,np.newaxis]*expdelta_kernel,0)

    
    ###############
    #Trigger Efficiency
    ###############
    #Function of calculated using true energies, before detector resolution effects
    if doEffs:
        n_Eee*=eff.trigEff(Ebin_ctr)

    ###############
    #Detector Resolution
    ###############
    if doDetRes:
        #Smear spectrum using resolution function
        gaus_kernel = gausInt_fast(Ebin_ctr[:,np.newaxis], sigma_ee(Ebin_ctr[:,np.newaxis],sigma0,B_1,A), Ebins)
        n_Eee = np.sum(n_Eee[:,np.newaxis]*gaus_kernel,0)

    
    ###############
    #Efficiencies
    ###############
    #Use the fitted, smoothed total cut efficiency curve
    if doEffs:
        n_Eee*=eff.cutEffFit(Ebin_ctr)

    return n_Eee*scale

###########################################################################
#Apply yield, resolution, efficiencies etc. to simulated composite NR hit data to get an average spectrum
##########################
#NOTE: This is for composite NR's like (n,gamma)
##########################
#
#Ebins: Bin edges in [eV]
#Evec: Array of hit Energy vectors
#dEvec: Array of hit deposited Energy vectors
#Yield: Either a constant or a Yield object with a .calc method as defined in R68_yield.py
#F: Fano factor
#scale: Scaling of spectrum 
#doDetRes: whether to apply detector resolution function
#fpeak: Vbias smearing effect parameter. Fraction of events which are NOT smeared. (1.0=no smearing) See EPot_Model.ipynb
#doEffs: Apply trigger and cut efficiency curves
#
#Returns: n_Eee: The binned spectra of events in units of [counts/bin]
#
def buildAvgSimSpectrum_ee_composite(Ebins, Evec, dEvec, Yield, F, scale, doDetRes=True, fpeak=1.0, doEffs=True,verbose=False):

    if verbose:
        print('In-function yield:', Yield.pars, 'In-function F_NR:', F)
    
    V=const.V
    G_NTL=const.G_NTL
    eps=const.eps
    F_ER=const.F
    
    #Params from Matt's Bkg resolution fit:
    #https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_panda:calibration#resolution_versus_energy
    sigma0=const.sigma0 #eV 
    B=const.B #This includes FANO
    #B_1=B-F*eps #Actually only true for V->inf (http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/Peak_Widths.pdf)
    B_1=B-F_ER*V**2/eps/(1+V/eps)**2
    A=const.A

    ###############
    #Binning
    ###############
    Ebin_ctr = (Ebins[:-1]+Ebins[1:])/2

    ###############
    #Yield, Fano
    ###############
    
    #Ionization energy available in each step
    if isinstance(Yield,(float,int)):
        Evec_eion = Evec*Yield - (Evec-dEvec)*Yield
    else:
        Evec_eion = Evec*Yield.calc(Evec) - (Evec-dEvec)*Yield.calc(Evec-dEvec)

    #Electron-equivalent
    EeeVec = (dEvec + Evec_eion*V/eps) / G_NTL
    #Fano resolution in eVee
    #TODO: This is not quite correct if dE/E<<1, but may be close enough since that's rare?
    SigmaEeeSqVec = Evec_eion*(V**2)*F/eps/(G_NTL**2)
    
    #energy==0 (no hit) may still generate Yield!=0, take care of that here
    SigmaEeeSqVec[Evec_eion<=0]=0.0
    
    #Sum hits assuming they're uncorrelated and reshape so we can broadcast against bin centers
    Eee=np.sum(EeeVec,1)[:,np.newaxis]
    SigmaEee = np.sqrt(np.sum(SigmaEeeSqVec,1))[:,np.newaxis]
    
    #Weighted spectra of Eee energies measured
    nvec_Eee = gausInt_fast(Eee, SigmaEee, Ebins)
    #Spectrum
    n_Eee = np.sum(nvec_Eee,0)

        
    ###############
    #Vbias Smearing
    ###############
    #Currently only for studying systematic uncertainties, smearing model is an over-estimation
    #See EPot_Model.ipynb
    if fpeak<1.0:
        #Smear spectrum using resolution function
        expdelta_kernel=expdelta_int(Ebin_ctr[:,np.newaxis],fpeak,4.20823831,Ebins)
        n_Eee = np.sum(n_Eee[:,np.newaxis]*expdelta_kernel,0)

    
    ###############
    #Trigger Efficiency
    ###############
    #Function of calculated using true energies, before detector resolution effects
    if doEffs:
        n_Eee*=eff.trigEff(Ebin_ctr)

    ###############
    #Detector Resolution
    ###############
    if doDetRes:
        #Smear spectrum using resolution function
        gaus_kernel = gausInt_fast(Ebin_ctr[:,np.newaxis], sigma_ee(Ebin_ctr[:,np.newaxis],sigma0,B_1,A), Ebins)
        n_Eee = np.sum(n_Eee[:,np.newaxis]*gaus_kernel,0)

    
    ###############
    #Efficiencies
    ###############
    #Analysis cut efficiencies
    #Use the fitted, smoothed total cut efficiency curve
    if doEffs:
        n_Eee*=eff.cutEffFit(Ebin_ctr)

    return n_Eee*scale


###########################################################################
#N e/h pairs function
#simplistic gaussian model for charge to get width and mean correct
#discretize a normal distribution, don't allow negative N
def getNeh(E,eps,F):
    #TMP
    #return np.random.normal((E/eps),np.sqrt(F*(E/eps)),np.shape(E))
    N=np.round(np.random.normal((E/eps),np.sqrt(F*(E/eps)),np.shape(E)))
    return np.maximum(0,N)

###########################################################################
#Select energies according to trigger efficiency
#Return whether energies passed the trigger
def passTrig(E):
    return np.random.random(E.shape)<eff.trigEff(E)

#Resolution functions

###########################################################################
#Detector Resolution function
def sigma_ee(E,sigma0,B,A):
    return np.sqrt(sigma0**2 + B*E + (A*E)**2)

###########################################################################
#Draw a value from the resolution distribution
#Option to include effect of OF bias at low energies
def getSmeared(E, F=0, doLowEbias=False, seed=None):
    #Params from Matt's Bkg resolution fit:
    #https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_panda:calibration#resolution_versus_energy
    sigma0=const.sigma0 #eV 
    B=const.B #This includes FANO
    # B_1=B-const.F*const.eps #Actually only true for V->inf (http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/Peak_Widths.pdf)
    B_1=B+(F-const.F)*const.V**2/const.eps/(1+const.V/const.eps)**2
    A=const.A
    
    #Ignore low energy bias for faster execution
    if not doLowEbias:
        if seed is not None:
            np.random.seed(seed)
        return np.random.normal(E,sigma_ee(E,sigma0,B_1,A),np.shape(E))
    else:
        #Low energy OF bias
        #From Nick's fit to PuBe noise wall
        N=3.83

        Esmeared=np.array([])
        for Ei in E:
            s_ee=sigma_ee(Ei,sigma0,B,A)
            vals=np.linspace(Ei-5*s_ee,Ei+5*s_ee,100)
            weights=P_OF_max(vals,Ei,s_ee,N)
            Esmeared=np.append(Esmeared,choices(vals, weights)[0])
        return Esmeared

###########################################################################
#OF resolution function
#Returns probability of Ahat given A and N indpendent bins when selecting the max Ahat
#See http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/OF_bias_theory.html
def P_OF_max(Ahat,A,sigma,N):
    term1=P_OF0(Ahat,A,sigma)*(0.5+0.5*erf(Ahat/np.sqrt(2)/sigma))**(N-1)
    term2A=(N-1)*P_OF0(Ahat,0,sigma)*(0.5+0.5*erf(Ahat/np.sqrt(2)/sigma))**(N-2)
    term2B=0.5+0.5*erf((Ahat-A)/np.sqrt(2)/sigma)
    return term1+term2A*term2B

###########################################################################
#PDF for OF0
def P_OF0(Ahat,A,sigma):
    return 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(Ahat-A)**2/2/sigma**2)

###########################################################################
#Integral of gaussian
def gausInt(mu, sigma, xlow, xhi):
    return np.where(sigma>0, 0.5*erf((mu-xlow)/np.sqrt(2)/sigma) - 0.5*erf((mu-xhi)/np.sqrt(2)/sigma), 1.0*((mu>xlow) & (mu<xhi)))

#Faster implementation of gausInt
#Uses lookup table for erf and smarter broadcasting steps
#Gives about 2x speed improvement
#mu: gaussian means. Should be float or shape == (N,1)
#sigma: gaussian widths. Should be float or shape == (N,1). Should all be sigma>=0
#bins: bin edges. Shape should be (M+1,)
def gausInt_fast(mu, sigma, bins):    
    #Precalculate these ratios
    #a=mu/np.sqrt(2)/sigma
    #b=bins/np.sqrt(2)/sigma
    
    #return np.where(sigma>0,0.5*erf_tab(a-b[:,:-1])-0.5*erf_tab(a-b[:,1:]), np.histogram(mu,bins)[0])
    #return np.where(sigma>0,0.5*erf_tab(a-b[:,:-1])-0.5*erf_tab(a-b[:,1:]), 1.0*((mu>bins[:-1]) & (mu<bins[1:])))
    
    if mu.ndim==0:
        
        if sigma==0:
            return 1.0*((mu>bins[:-1]) & (mu<bins[1:]))
        else:
            erf_res=erf_tab((mu-bins)/np.sqrt(2)/sigma)
            return -0.5*np.diff(erf_res)
        
    elif mu.ndim==2:
        ###Avoid divide by 0
        #Preallocate the output array
        #Assumes shapes described above
        N=mu.shape[0]
        M=bins.shape[0]-1
        result=np.zeros((N,M))

        #Where is sigma zero?
        sig0=(np.argwhere(sigma.reshape((-1))==0)).reshape((-1))
        signz=(np.argwhere(sigma.reshape((-1))!=0)).reshape((-1))

        #Handle sigma==0
        result[sig0]=1.0*((mu[sig0]>bins[:-1]) & (mu[sig0]<bins[1:]))

        #Handle sigma>0
        #Shouldn't be any sigma<0...
        erff=erf_tab(mu[signz]/np.sqrt(2)/sigma[signz]-bins/np.sqrt(2)/sigma[signz])

        result[signz]=-0.5*np.diff(erff)

        return result
    else:
        print('Error in R68_spec_tools.gausInt_fast. Mu has wrong dimensionality: ',mu.ndim)
        return None

    
#Lookup table for erf
x_erf = np.linspace(-3, 3, 1001)    #There are better sampling points, but this is sufficient to ~0.001%
y_erf = [erf(_) for _ in x_erf ]
def erf_tab(z):
    return np.interp(z, x_erf, y_erf, left=-1,right=1)
    
###########################################################################
#Smearing function for fringing field effect

#Integral of exponential PDF:P(x) = A*exp(B*(x-x0)/x0) + fpeak*delta(x-x0)
#A is normalized such that the total integral is 1
#
#bins: bin edges to integrate between
#x0: position of peak which we want to smear down
#fpeak: fraction of events in the delta function peak, from fringing field simulations
#B: fall rate [unitless], from fringing field simulations
def expdelta_int(x0,fpeak,B, bins):
    
    arg=B*(bins-x0)/x0

    #Integration bounds are (0,x0)
    arg[bins>x0]=0
    arg[bins<np.zeros_like(x0)]=-B
    
    #Calulate just once for speed
    below_x0=bins<x0
    
    return ((1-fpeak)/(1-np.exp(-B)))*(np.exp(arg[:,1:]) - np.exp(arg[:,:-1])) + (below_x0[:,:-1] & ~below_x0[:,1:])*fpeak

###########################################################################
#Plot of individual and combined simulated spectra and observed spectrum
#If N_nr, N_er, N_ng are more than 1-d, then the first row is assumed to be 
#  the best-fit spectrum, the second is the upper 1-sigma line and the last is the lower 1-sigma line.
#  Then we plot the best fit and shade between the upper and lower lines.
# For backwards compatibility with some notebooks, if N_tot is not included, it, 
#  and the corresponding uncertainty, are calculated from the sum of N_nr, N_er, and N_ng
def plotSpectra(E_bins, N_nr, N_er, N_ng, N_meas, dN_meas, N_tot=None, xrange=(0,1e3), yrange=(0,1e-2), yscale='linear', thresh=None, axis=None, wLeg=True, grid=True, wResidual=False, yrange_res=(-1e-3,1e-3)):
    
    #######################
    #Plotting styles
    import os
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
    #######################
    
    if axis is None:
        fig_w=9
        if wResidual:
            #fig,axis = plt.subplots(1,2,figsize=(fig_w, 2*fig_w*(.75)), sharex=True)
            #TODO Make Residual narrow
            fig = plt.figure(figsize=(fig_w, 1.5*fig_w*(.75)))
            gs = fig.add_gridspec(16, 1)
            ax0 = fig.add_subplot(gs[:10, :])
            ax1 = fig.add_subplot(gs[12:, :])
            axis=[ax0,ax1]
        else:
            fig,axis = plt.subplots(1,1,figsize=(fig_w, fig_w*(.75)))
    else:
        if wResidual:
            print()#TODO Make sure axis is an array with 2 axis entries

    assert N_nr.shape == N_er.shape == N_ng.shape
    if N_nr.ndim==2 and N_nr.shape[0]==3:
        dolimits=True
        #We've got a median and 2 limit spectra
        N_nr_best=N_nr[0]
        N_er_best=N_er[0]
        N_ng_best=N_ng[0]
        
        N_nr_up=N_nr[1]
        N_er_up=N_er[1]
        N_ng_up=N_ng[1]
        
        N_nr_down=N_nr[2]
        N_er_down=N_er[2]
        N_ng_down=N_ng[2]
        
        if N_tot is None:
            N_tot_best=N_nr_best+N_er_best+N_ng_best
            #TODO: What's the right way to combine these uncertainties?
            #Hard to say at this point if we don't know how they're correlated
            N_tot_up=N_nr_up+N_er_up+N_ng_up
            N_tot_down=N_nr_down+N_er_down+N_ng_down
        else:
            N_tot_best=N_tot[0]
            N_tot_up=N_tot[1]
            N_tot_down=N_tot[2]
        
    else:
        #Assume it's just a single best fit spectrum
        dolimits=False
        N_nr_best=N_nr
        N_er_best=N_er
        N_ng_best=N_ng
        
        if N_tot is None:
            N_tot_best=N_nr_best+N_er_best+N_ng_best
        else:
            N_tot_best=N_tot
        
    #Plot Fit
    dE = E_bins[1]-E_bins[0]
    Ec = (E_bins[:-1] + E_bins[1:]) / 2
    
    if wResidual:
        ax=axis[0]
    else:
        ax=axis
        
    ax.step(Ec,N_nr_best, where='mid',color='k', linestyle='-', \
             label='NR, no Capture', linewidth=2)
    
    ax.step(Ec,N_er_best, where='mid',color='r', linestyle='-', \
             label='ER, no Capture', linewidth=2)

    ax.step(Ec,N_ng_best, where='mid',color='b', linestyle='-', \
             label='(n,gamma)', linewidth=2)

    ax.step(Ec,N_tot_best, where='mid',color='g', linestyle='-',\
            label='All Sims', linewidth=2)


    if dolimits:
        ax.fill_between(Ec,N_nr_down,N_nr_up, step='mid', color='k', alpha=0.5)
        ax.fill_between(Ec,N_er_down,N_er_up, step='mid', color='r', alpha=0.5)
        ax.fill_between(Ec,N_ng_down,N_ng_up, step='mid', color='b', alpha=0.5)
        ax.fill_between(Ec,N_tot_down,N_tot_up, step='mid', color='g', alpha=0.5)
    
    #Plot Data
    if np.asarray(dN_meas).ndim==2:
        #asymmetric
        dN_meas_hi=dN_meas[0]
        dN_meas_low=dN_meas[1]
    else:
        #symmetric
        dN_meas_hi=dN_meas
        dN_meas_low=dN_meas
    ax.errorbar(Ec,N_meas,yerr=[dN_meas_hi,dN_meas_low], marker='o', markersize=6, \
                 ecolor='k',color='k', linestyle='none', label='Data-Bkg', linewidth=2)
    
    if thresh is not None:
        ax.axvline(thresh, color='m', linestyle='--', linewidth=2, label='Threshold')
    
    #Plot Residual
    if wResidual:
        axis[1].step(Ec,np.zeros_like(Ec), where='mid',color='g', linestyle='-', linewidth=2)
        axis[1].errorbar(Ec,N_meas-N_tot_best,yerr=[dN_meas_hi,dN_meas_low], marker='o', markersize=6, \
                 ecolor='k',color='k', linestyle='none', linewidth=2)
        axis[1].set_ylim(*yrange_res)
        axis[1].set_xlim(*xrange)
    
    ax.set_yscale(yscale)
    ax.set_xlim(*xrange)
    ax.set_ylim(*yrange)
    #ax.set_xlabel('total deposited energy [eV$_{\\mathrm{ee}}$]',**axis_font)
    ax.set_xlabel('Energy [eV$_{\\mathrm{ee}}$]',**axis_font)
    #axis.set_ylabel('counts',**axis_font)
    #ax.set_ylabel('Events/bin/s',**axis_font)
    ax.set_ylabel('Rate [1/bin/s]',**axis_font)
    if grid:
        ax.grid(True)
        ax.yaxis.grid(True,which='minor',linestyle='--')
    if wLeg:
        ax.legend(loc=1,prop={'size':22})

    for side in ['top','bottom','left','right']:
      ax.spines[side].set_linewidth(2)

    plt.tight_layout()
    
    return axis