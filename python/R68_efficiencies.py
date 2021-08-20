#R68 Efficiency functions
#All energy units are eVee
import numpy as np

#Trigger efficiency
###################
#PuBe
#Note, this is a function of true pulse energy, before OF resolution effects
with open('data/r68_trigger_eff_1keV.txt') as fefftrig:
    data_efftrig = np.asarray([x.split() for x in fefftrig.readlines()[1:]],dtype=np.float)

def trigEff(E):
    #https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_trigger
    #v3
    #return erf(E/(np.sqrt(2.)*30.86))
    #v4
    #return 1-np.exp(-((E/57.3015)**1.07219))

    #Final
    return np.interp(E,data_efftrig[:,0],data_efftrig[:,1])
    
#Return uncertainties on the trigger efficiency.
#These come from the stdev across multiple trigger sims which dominated statistical uncertainties in a given sim
def dtrigEff(E):
    return np.interp(E,data_efftrig[:,0],data_efftrig[:,2])


#Write efficiencies
###################
#Note these are NOT included in the livetimes calculated in R68_load
eff_write = 0.617 #https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_rateandlivetime#iii_read_efficiency
deff_write = 0.004
eff_write_bkg = 0.815
deff_write_bkg = 0.004

#Cut Efficiencies
#################
eff_tail = 0.8197 #http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/Run68_CutEff_1/
deff_tail = 0.0013
eff_pileup = 0.9651 #http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/Run68_CutEff_2/
deff_pileup = 0.0013
#Replaced with energy-dependent version
#eff_trigburst = 0.9887 #http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/Run68_RateCut_pt2/
#deff_trigburst = 0.0013 #Estimated from the out-of band passage which has uncertainty ~sqrt(2/N) and N~1.2e6

eff_tail_bkg = 0.8875
deff_tail_bkg = 0.0026
eff_pileup_bkg = 0.9779
deff_pileup_bkg = 0.0040
#eff_trigburst_bkg = 0.9992
#deff_trigburst_bkg = 0.0013

#http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/Run68_RateCut_pt2/
#Replaced with dN cut: burst_cut.ipynb
def trigburstEff(E):
    #eff=np.zeros_like(E)
    #eff[E<=37.01]=0.9863
    #eff[(E>37.01) & (E<=148.9)]=0.9709
    #eff[(E>148.9) & (E<=376.2)]=0.9825
    #eff[E>376.2]=0.9863
    #return eff
    eff=np.ones_like(E)
    eff[(E>=50)&(E<1000)]=0.8933
    return eff
    
def dtrigburstEff(E):
    #return 0.0021*np.ones_like(E)
    deff=np.zeros_like(E)
    deff[(E>=50)&(E<1000)]=0.0006
    return deff

#No dN cut needed for bkg
#def trigburstEff_bkg(E):
#    return 0.9992*np.ones_like(E)
#def dtrigburstEff_bkg(E):
#    return 0.0021*np.ones_like(E)

#Estimated maximum leakage spectrum of burst events
# Fit performed in burst_cut.ipynb
# Assumes uniform bins and full measured PuBe stats
def dtrigburstLeak(E):
    def fburst(E,A,E0,Ec,Efall):
        hump=(E>E0)*(np.exp(-(E-E0)/Efall)-np.exp(-(E-E0)/Ec))
        hump=A*hump/(np.max(hump))
        return hump
    dE0=10 #[eV] Energy binning in the fit
    dE=E[1]-E[0] #[eV] Scale by binning ratio to maintain (roughly) the same integral
    return (dE/dE0)*fburst(E, 98.74073987, 51.7184786, 3.09597772, 39.05670093)

#Spikey Cut Efficiency
with open('data/r68_PuBe_cspike_eff_1keV.txt') as feffspike:
    data_effspike = np.asarray([x.split() for x in feffspike.readlines()[1:]],dtype=np.float)

with open('data/eff_spikeSim_r68_bkg.csv') as feffspike_bkg:
    data_effspike_bkg = np.asarray([x.split(',') for x in feffspike_bkg.readlines()],dtype=np.float)
    
def spikeEff(E):
    #http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/Run68_SpikeEff/
    return np.interp(E,data_effspike[:,0],data_effspike[:,1])

#Return upper and lower spike cut uncertainties
def dspikeEff(E):
    dup = np.interp(E,data_effspike[:,0],data_effspike[:,2])
    dlow = np.interp(E,data_effspike[:,0],data_effspike[:,3])
    return np.stack((dup,dlow))

def spikeEff_bkg(E):
    return np.interp(E,data_effspike_bkg[:,0],data_effspike_bkg[:,1])
def dspikeEff_bkg(E):
    dup = np.interp(E,data_effspike_bkg[:,0],data_effspike_bkg[:,2])
    dlow = dup
    return np.stack((dup,dlow))

#Chisq Cut Efficiency
with open('data/r68_PuBe_cchit_eff_1keV.txt') as feffchi:
    data_effchi = np.asarray([x.split() for x in feffchi.readlines()[1:]],dtype=np.float)

with open('data/r68_bkg_cchit_eff_1keV.txt') as feffchi_bkg:
    data_effchi_bkg = np.asarray([x.split() for x in feffchi_bkg.readlines()[1:]],dtype=np.float)

def chisqEff(E):
    #v1
    #eff_chisq = 0.645 #http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/Run68_CutEff_2/
    
    #v2 w/Energy dependence
    #http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/Run68_CutEff_2/
    return np.interp(E,data_effchi[:,0],data_effchi[:,1])

#Return upper and lower chisq cut uncertainties
def dchisqEff(E):
    dup = np.interp(E,data_effchi[:,0],data_effchi[:,2])
    dlow = np.interp(E,data_effchi[:,0],data_effchi[:,3])
    return np.stack((dup,dlow))

def chisqEff_bkg(E):
    return np.interp(E,data_effchi_bkg[:,0],data_effchi_bkg[:,1])
def dchisqEff_bkg(E):
    dup = np.interp(E,data_effchi_bkg[:,0],data_effchi_bkg[:,2])
    dlow = np.interp(E,data_effchi_bkg[:,0],data_effchi_bkg[:,3])
    return np.stack((dup,dlow))

#Return the total CUT efficiency curve (i.e. no trigger or write efficiency)
def cutEff(E):
    return eff_tail*eff_pileup*trigburstEff(E)*spikeEff(E)*chisqEff(E)
def cutEff_bkg(E):
    #return eff_tail_bkg*eff_pileup_bkg*trigburstEff_bkg(E)*spikeEff_bkg(E)*chisqEff_bkg(E)
    return eff_tail_bkg*eff_pileup_bkg*spikeEff_bkg(E)*chisqEff_bkg(E)


#Return the upper and lower total cut uncertainties
#Adding asymm errors in quadrature is apparently naughty (https://www.slac.stanford.edu/econf/C030908/papers/WEMT002.pdf)
#But they're not that asymmetric, so it's probably fine
def dcutEff(E):
    
    dupsq = (deff_write/eff_write)**2 + (deff_tail/eff_tail)**2 + (deff_pileup/eff_pileup)**2 + \
    (dtrigburstEff(E)/trigburstEff(E))**2 + (dspikeEff(E)[0]/spikeEff(E))**2 + (dchisqEff(E)[0]/chisqEff(E))**2
    #dup = np.sqrt(dupsq/(cutEff(E)**2)) This isn't right at all!
    dup = cutEff(E)*np.sqrt(dupsq)
    
    dlowsq = (deff_write/eff_write)**2 + (deff_tail/eff_tail)**2 + (deff_pileup/eff_pileup)**2 + \
    (dtrigburstEff(E)/trigburstEff(E))**2 + (dspikeEff(E)[1]/spikeEff(E))**2 + (dchisqEff(E)[1]/chisqEff(E))**2
    #dlow = np.sqrt(dlowsq/(cutEff(E)**2))
    dlow = cutEff(E)*np.sqrt(dlowsq)
    
    return np.stack((dup,dlow))

def dcutEff_bkg(E):
    
    dupsq = (deff_write_bkg/eff_write_bkg)**2 + (deff_tail_bkg/eff_tail_bkg)**2 + (deff_pileup_bkg/eff_pileup_bkg)**2 + \
    + (dspikeEff_bkg(E)[0]/spikeEff_bkg(E))**2 + (dchisqEff_bkg(E)[0]/chisqEff_bkg(E))**2
    #dup = np.sqrt(dupsq/(cutEff_bkg(E)**2))
    dup = cutEff_bkg(E)*np.sqrt(dupsq)
    
    dlowsq = (deff_write_bkg/eff_write_bkg)**2 + (deff_tail_bkg/eff_tail_bkg)**2 + (deff_pileup_bkg/eff_pileup_bkg)**2 + \
    + (dspikeEff_bkg(E)[1]/spikeEff_bkg(E))**2 + (dchisqEff_bkg(E)[1]/chisqEff_bkg(E))**2
    #dlow = np.sqrt(dlowsq/(cutEff_bkg(E)**2))
    dlow = cutEff_bkg(E)*np.sqrt(dlowsq)
    
    return np.stack((dup,dlow))

#Here is the total efficiency including the trigger efficiency
def allEff(E):
    return eff_tail*eff_pileup*trigburstEff(E)*spikeEff(E)*chisqEff(E)*trigEff(E)
def allEff_bkg(E):
    #return eff_tail_bkg*eff_pileup_bkg*trigburstEff_bkg(E)*spikeEff_bkg(E)*chisqEff_bkg(E)
    return eff_tail_bkg*eff_pileup_bkg*spikeEff_bkg(E)*chisqEff_bkg(E)*trigEff(E)


#Return the upper and lower total cut uncertainties
#Adding asymm errors in quadrature is apparently naughty (https://www.slac.stanford.edu/econf/C030908/papers/WEMT002.pdf)
#But they're not that asymmetric, so it's probably fine
def dallEff(E):
    
    dupsq = (deff_write/eff_write)**2 + (deff_tail/eff_tail)**2 + (deff_pileup/eff_pileup)**2 + \
    (dtrigburstEff(E)/trigburstEff(E))**2 + (dspikeEff(E)[0]/spikeEff(E))**2 + (dchisqEff(E)[0]/chisqEff(E))**2 + \
    (dtrigEff(E)[0]/trigEff(E))**2
    #dup = np.sqrt(dupsq/(cutEff(E)**2)) This isn't right at all!
    dup = cutEff(E)*np.sqrt(dupsq)
    
    dlowsq = (deff_write/eff_write)**2 + (deff_tail/eff_tail)**2 + (deff_pileup/eff_pileup)**2 + \
    (dtrigburstEff(E)/trigburstEff(E))**2 + (dspikeEff(E)[1]/spikeEff(E))**2 + (dchisqEff(E)[1]/chisqEff(E))**2 + \
    (dtrigEff(E)[1]/trigEff(E))**2
    #dlow = np.sqrt(dlowsq/(cutEff(E)**2))
    dlow = cutEff(E)*np.sqrt(dlowsq)
    
    return np.stack((dup,dlow))

def dallEff_bkg(E):
    
    dupsq = (deff_write_bkg/eff_write_bkg)**2 + (deff_tail_bkg/eff_tail_bkg)**2 + (deff_pileup_bkg/eff_pileup_bkg)**2 + \
    + (dspikeEff_bkg(E)[0]/spikeEff_bkg(E))**2 + (dchisqEff_bkg(E)[0]/chisqEff_bkg(E))**2
    #dup = np.sqrt(dupsq/(cutEff_bkg(E)**2))
    dup = cutEff_bkg(E)*np.sqrt(dupsq)
    
    dlowsq = (deff_write_bkg/eff_write_bkg)**2 + (deff_tail_bkg/eff_tail_bkg)**2 + (deff_pileup_bkg/eff_pileup_bkg)**2 + \
    + (dspikeEff_bkg(E)[1]/spikeEff_bkg(E))**2 + (dchisqEff_bkg(E)[1]/chisqEff_bkg(E))**2
    #dlow = np.sqrt(dlowsq/(cutEff_bkg(E)**2))
    dlow = cutEff_bkg(E)*np.sqrt(dlowsq)
    
    return np.stack((dup,dlow))


#Here is a smoothed efficiency function to fit to, and use in place of, the observed efficiency curves.
def effFit_func(x,x0,sigma,a,b):
    return (x>x0)*(a+b*x)*(1-np.exp(-(x-x0)/sigma))

#Best fit values to total cut efficiency as calculated in R68_eff_plot.ipynb
#PuBe
#def cutEffFit(E):
#    return effFit_func(E, 2.51808946e+01, 1.49032295e+01, 4.66827422e-01, 4.43609440e-05)
#def dcutEffFit(E):
#    return 0.03430284917913029

def cutEffFit(E):
    eff=np.ones_like(E)
    eff[(E>=50)&(E<1000)]=effFit_func(E[(E>=50)&(E<1000)], 2.16419228e+01, 1.61938240e+01, 4.28262762e-01, 3.17895378e-05)
    eff[E>=1000]=0.51633498
    return eff
def dcutEffFit(E):
    deff=np.zeros_like(E)
    deff[(E>=50)&(E<1000)]=0.02975340974950752
    deff[E>=1000]=0.01294113
    return deff

#Bkg
#def cutEffFit_bkg(E):
#    return effFit_func(E, 3.62217815e+01, 1.19451268e+01, 7.13281258e-01, 7.41533280e-06)
#def dcutEffFit_bkg(E):
#    return 0.16433241774395257
def cutEffFit_bkg(E):
    return effFit_func(E, 4.23508507e+01,  7.61590223e+00,  7.14511426e-01, -5.10367727e-06)
def dcutEffFit_bkg(E):
    return 0.16643841685131444


#Best fit values to total cut efficiency as calculated in R68_all_eff_plot.ipynb (including trigger eff)
#PuBe
def allEffFit(E):
    eff=np.ones_like(E)
    eff[(E>=50)&(E<1000)]=effFit_func(E[(E>=50)&(E<1000)], 1.99011656e+01, 3.39505540e+01, 4.07808092e-01, 6.02809942e-05)
    eff[E>=1000]=0.51633498
    return eff
def dallEffFit(E):
    deff=np.zeros_like(E)
    deff[(E>=50)&(E<1000)]=0.03230678598469722
    deff[E>=1000]=0.01294113
    return deff

#Bkg
def allEffFit_bkg(E):
    return effFit_func(E, 2.12459234e+01,  3.66441228e+01,  7.02002489e-01, 1.40780632e-06)
def dallEffFit_bkg(E):
    return 0.15414768212039148 
