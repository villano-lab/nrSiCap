###########################################################################
#Methods for working with spectra in R68
###########################################################################
import constants as const
import R68_efficiencies as eff
from scipy.special import erf
exec(open("nb_setup.py").read()) #Is there a better way to do this?

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
#
#Returns: (n_Eee_nr, n_Eee_er, n_Eee_ng) The binned spectra of NR, ER, and (n,gamma) events in units of [counts/bin]

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
def buildAvgSimSpectra_ee(Ebins, Evec_nr, Evec_er, Evec_ng, dEvec_ng, Yield, F_NR, scale_g4=1, scale_ng=1, doDetRes=True, seed=1):

    V=const.V
    G_NTL=const.G_NTL
    eps=const.eps
    F=const.F
    
    #Params from Matt's Bkg resolution fit:
    #https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_panda:calibration#resolution_versus_energy
    sigma0=const.sigma0 #eV 
    B=const.B #This includes FANO
    B_1=B-F*eps
    A=const.A
    
    np.random.seed(seed)

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
    nvec_Eee_nr = gausInt(Eee_nr, SigmaEee_nr, Ebins[:-1], Ebins[1:])
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
    nvec_Eee_er = gausInt(Eee_er, SigmaEee_er, Ebins[:-1], Ebins[1:])
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
    nvec_Eee_ng = gausInt(Eee_ng, SigmaEee_ng, Ebins[:-1], Ebins[1:])
    #(n,gamma) spectrum
    n_Eee_ng = np.sum(nvec_Eee_ng,0)

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
        gaus_kernel = gausInt(Ebin_ctr[:,np.newaxis], sigma_ee(Ebin_ctr[:,np.newaxis],sigma0,B_1,A), Ebins[:-1], Ebins[1:])
        
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
def getSmeared(E, doLowEbias=False):
    #Params from Matt's Bkg resolution fit:
    #https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_panda:calibration#resolution_versus_energy
    sigma0=const.sigma0 #eV 
    B=const.B #This includes FANO
    B_1=B-const.F*const.eps
    A=const.A
    
    #Ignore low energy bias for faster execution
    if not doLowEbias:
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
    return np.where(sigma>0, 0.5*erf((mu-xlow)/np.sqrt(2*sigma**2)) - 0.5*erf((mu-xhi)/np.sqrt(2*sigma**2)), 1.0*((mu>xlow) & (mu<xhi)))

###########################################################################
#Plot of individual and combined simulated spectra and observed spectrum
def plotSpectra(E_bins, N_nr, N_er, N_ng, N_meas, dN_meas, xrange=(0,1e3), yrange=(0,1e-2), yscale='linear', thresh=None):
    
    fig_w=9
    fig,axes = plt.subplots(1,1,figsize=(fig_w, fig_w*(.75)), sharex=True)
    ax1 = axes

    dE = E_bins[1]-E_bins[0]
    
    Ec = (E_bins[:-1] + E_bins[1:]) / 2
    
    ax1.step(Ec,N_nr, where='mid',color='k', linestyle='-', \
             label='NR, no Capture', linewidth=2)

    ax1.step(Ec,N_er, where='mid',color='r', linestyle='-', \
             label='ER, no Capture', linewidth=2)

    ax1.step(Ec,N_ng, where='mid',color='b', linestyle='-', \
             label='(n,gamma)', linewidth=2)

    ax1.step(Ec,N_nr+N_er+N_ng, where='mid',color='g', linestyle='-', \
             label='All Sims', linewidth=2)

    ax1.errorbar(Ec,N_meas,yerr=[dN_meas,dN_meas], marker='o', markersize=6, \
                 ecolor='k',color='k', linestyle='none', label='Data-Bkg', linewidth=2)
    
    if thresh is not None:
        ax1.axvline(thresh, color='m', linestyle='--', linewidth=2, label='Threshold')
    
    ax1.set_yscale(yscale)
    ax1.set_xlim(*xrange)
    ax1.set_ylim(*yrange)
    ax1.set_xlabel('total deposited energy [eV$_{\\mathrm{ee}}$]',**axis_font)
    #ax1.set_ylabel('counts',**axis_font)
    ax1.set_ylabel('Events/bin/s',**axis_font)
    ax1.grid(True)
    ax1.yaxis.grid(True,which='minor',linestyle='--')
    ax1.legend(loc=1,prop={'size':22})

    for axis in ['top','bottom','left','right']:
      ax1.spines[axis].set_linewidth(2)

    plt.tight_layout()