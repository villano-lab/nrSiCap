#R68 Efficiency functions
import numpy as np

#Trigger efficiency
#Note, this is a function of true pulse energy, before OF resolution effects
fefftrig = open('data/r68_trigger_eff_1keV.txt')
data_efftrig = np.asarray([x.split() for x in fefftrig.readlines()[1:]],dtype=np.float)
fefftrig.close()

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


#Cut Efficiencies
eff_write = 0.617 #https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_rateandlivetime#iii_read_efficiency
deff_write = 0.004
eff_tail = 0.8197 #http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/Run68_CutEff_1/
deff_tail = 0.0013
eff_pileup = 0.9651 #http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/Run68_CutEff_2/
deff_pileup = 0.0013
eff_trigburst = 0.9887 #http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/Run68_RateCut_pt2/
deff_trigburst = 0.0013 #Estimated from the out-of band passage which has uncertainty ~sqrt(2/N) and N~1.2e6

#Spikey Cut Efficiency
feffspike = open('data/r68_cspike_eff_1keV.txt')
data_effspike = np.asarray([x.split() for x in feffspike.readlines()[1:]],dtype=np.float)
feffspike.close()

def spikeEff(E):
    #http://www.hep.umn.edu/cdms/cdms_restricted/K100/analysis/Run68_SpikeEff/
    return np.interp(E,data_effspike[:,0],data_effspike[:,1])

#Return upper and lower spike cut uncertainties
def dspikeEff(E):
    dup = np.interp(E,data_effspike[:,0],data_effspike[:,2])
    dlow = np.interp(E,data_effspike[:,0],data_effspike[:,3])
    return np.stack((dup,dlow))

#Chisq Cut Efficiency
feffchi = open('data/r68_cchit_eff_1keV.txt')
data_effchi = np.asarray([x.split() for x in feffchi.readlines()[1:]],dtype=np.float)
feffchi.close()

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

#Return the total cut efficiency curve
def cutEff(E):
    return eff_write*eff_tail*eff_pileup*eff_trigburst*spikeEff(E)*chisqEff(E)

#Return the upper and lower total cut uncertainties
#Adding asymm errors in quadrature is apparently naughty (https://www.slac.stanford.edu/econf/C030908/papers/WEMT002.pdf)
#But they're not that asymmetric, so it's probably fine
def dcutEff(E):
    
    dupsq = (deff_write/eff_write)**2 + (deff_tail/eff_tail)**2 + (deff_pileup/eff_pileup)**2 + \
    (deff_trigburst/eff_trigburst)**2 + (dspikeEff(E)[0]/spikeEff(E))**2 + (dchisqEff(E)[0]/chisqEff(E))**2
    dup = np.sqrt(dupsq/(cutEff(E)**2))
    
    dlowsq = (deff_write/eff_write)**2 + (deff_tail/eff_tail)**2 + (deff_pileup/eff_pileup)**2 + \
    (deff_trigburst/eff_trigburst)**2 + (dspikeEff(E)[1]/spikeEff(E))**2 + (dchisqEff(E)[1]/chisqEff(E))**2
    dlow = np.sqrt(dlowsq/(cutEff(E)**2))
    
    return np.stack((dup,dlow))
