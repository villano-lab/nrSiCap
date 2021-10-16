#Functions to load measured and simulated data for R68
import numpy as np

############################################################################
#Load measured R68 PuBe and Bkg data
#keVmax={2,10}
def load_measured(keVmax=2,verbose=False):
    if verbose:
        print('Loading Measured Data...')
    #Only low energy rate cut
    #with open('data/r68_n125_PuBe_cgood_final_PTOFkeV_2keV_scan_fmt.txt') as fpube:
    #Extended energy rate cut, nominal burst cut values
    #with open(f'data/r68_n125_PuBe_cgood_cRate{cLERate}_PTOFkeV_2keV_scan_fmt.txt') as fpube:
    
    #improved dN cut for bursts
    with open(f'data/r68_n125_PuBe_cok_cdN_PTOFkeV_{keVmax}keV_scan_fmt.txt') as fpube:
        dE = np.asarray([x.split() for x in fpube.readlines()],dtype=float)
        
    #with open('data/r68_n125_bkg_cgood_final_PTOFkeV_2keV_scan_fmt.txt') as fbknd:
    #No burst cut
    with open(f'data/r68_n125_bkg_cok_PTOFkeV_{keVmax}keV_scan_fmt.txt') as fbknd:
        dbE = np.asarray([x.split() for x in fbknd.readlines()],dtype=float)

    #Measured event energies in [eV] for PuBe (dE) and background (dbE)
    E_PuBe = dE[:,1]*1e3
    E_Bkg = dbE[:,1]*1e3
    if verbose:
        print('PuBe events: ', np.shape(E_PuBe))
        print('Bkg events: ', np.shape(E_Bkg))

    #Measured data live time estimates
    #https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_rateandlivetime#iii_read_efficiency
    tlive_PuBe = 97.9*3600 #[s] This is the raw livetime, without write efficiency effects
    
    #tlive_bkg = tlive_PuBe*193/640 # Naive scaling by number of series
    #tlive_bkg = 24.1*3600 #[s]
    tlive_bkg = 29.5*3600 #[s] Use the raw livetime, without applying write efficiency effects
    
    return {"Bkg":{"E":E_Bkg, "tlive":tlive_bkg},"PuBe":{"E":E_PuBe, "tlive":tlive_PuBe}}

############################################################################
#load simulated Geant4 data
#load_frac: fraction of simulated events to actually load
def load_G4(load_frac=1.0,verbose=False):
    if verbose:
        print('Loading Geant4 Data...')
    #===============to suppress h5py warning see:
    #https://github.com/h5py/h5py/issues/961
    import warnings
    warnings.simplefilter(action='ignore', category=FutureWarning)
    import h5py
    warnings.resetwarnings()
    import pandas as pd

    #f_nr_nocap = h5py.File("/data/chocula/villaa/k100Sim_Data/captureCalhdf5/R68_gdirect_testskim_superhighstat_nocap.h5","r")
    f_nr_nocap = h5py.File("data/R68_gdirect_testskim_stupidhighstat_nocap_nr.h5","r")
    data_nr_nocap = f_nr_nocap['geant4/hits']

    #f_er_nocap = h5py.File("/data/chocula/villaa/k100Sim_Data/captureCalhdf5/R68_gdirect_testskim_superhighstat_nocap_er_lowe.h5","r")
    f_er_nocap = h5py.File("data/R68_gdirect_testskim_stupidhighstat_nocap_er_lowe.h5","r")
    data_er_nocap = f_er_nocap['geant4/hits']
    if verbose:
        print(np.shape(data_nr_nocap))
        print(np.shape(data_er_nocap))

    #https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_n125:full_signal_fit&#efficiencies_and_live_time
    #tlive_g4 = 18.9*3600*load_frac #s
    tlive_g4 = 48.0*3600*load_frac #s

    #now make a dataframe with the restricted data
    #Columns are:
    #cols=['EV', 'DT', 'TS', 'P', 'Type', 'E1', 'D3', 'PX3', 'PY3', 'PZ3', 'X3', 'Y3', 'Z3',
    #      'time3', 'PX1', 'PY1', 'PZ1', 'X1', 'Y1', 'Z1', 'time1', 'nCap']

    cols=['EV', 'D3', 'X3', 'Y3', 'Z3','time3', 'nCap']

    sel_names=['EV', 'D3', 'X3', 'Y3', 'Z3', 'time3', 'nCap']#Select these variables
    sel=[cols.index(i) for i in sel_names]

    nr_nocap_data = data_nr_nocap[:,sel]
    nr_nocap_dataframe = pd.DataFrame(data=nr_nocap_data, columns=sel_names)
    er_nocap_data = data_er_nocap[:,sel]
    er_nocap_dataframe = pd.DataFrame(data=er_nocap_data, columns=sel_names)

    #need unique event numbers in case of (rare) duplicate 'EV's
    nr_nocap_evnew=np.cumsum(np.diff(nr_nocap_data[:,0],prepend=nr_nocap_data[0,0]).astype(bool).astype(float))
    nr_nocap_dataframe.insert(0,'EVnew',nr_nocap_evnew)

    er_nocap_evnew=np.cumsum(np.diff(er_nocap_data[:,0],prepend=er_nocap_data[0,0]).astype(bool).astype(float))
    er_nocap_dataframe.insert(0,'EVnew',er_nocap_evnew)

    ##Group hits into vectors for each thrown event
    #These loops can take a while
    import time

    groupbyvec=['EVnew']

    #NR, no capture
    if verbose:
        print("Loading NRs...")
    start = time.time()
    nev_nr_nocap = int(nr_nocap_dataframe.EVnew.nunique()*load_frac)
    max_vec_nr_nocap = np.max(nr_nocap_dataframe.groupby(groupbyvec).size()[:nev_nr_nocap-1])

    evec_nr_nocap = np.zeros((nev_nr_nocap,max_vec_nr_nocap))#Hit energies
    nhit_nr_nocap = np.zeros((nev_nr_nocap,1))#Number of hits

    nr_nocap_grouped=nr_nocap_dataframe.groupby(groupbyvec).agg({'D3':list})[:nev_nr_nocap-1]
    for iev,d3 in enumerate(nr_nocap_grouped.D3):
        evec_nr_nocap[iev,0:len(d3)] = d3
        nhit_nr_nocap[iev] = len(d3)
    evec_nr_nocap*=1e6 #[MeV]->[eV]

    end = time.time()
    if verbose:
        print(round((end - start)/60.,1),' min')

    #ER, no capture
    if verbose:
        print("Loading ERs...")
    start = time.time()
    nev_er_nocap = int(er_nocap_dataframe.EVnew.nunique()*load_frac)
    max_vec_er_nocap = np.max(er_nocap_dataframe.groupby(groupbyvec).size()[:nev_er_nocap-1])

    evec_er_nocap = np.zeros((nev_er_nocap,max_vec_er_nocap))#Hit energies
    nhit_er_nocap = np.zeros((nev_er_nocap,1))#Number of hits

    er_nocap_grouped=er_nocap_dataframe.groupby(groupbyvec).agg({'D3':list})[:nev_er_nocap-1]
    for iev,d3 in enumerate(er_nocap_grouped.D3):
        evec_er_nocap[iev,0:len(d3)] = d3 #[eV]
        nhit_er_nocap[iev] = len(d3)
    evec_er_nocap*=1e6 #[MeV]->[eV]
    end = time.time()
    if verbose:
        print(round((end - start)/60.,1),' min')
    
    return {"ER":{"E":evec_er_nocap, "N":nhit_er_nocap, "tlive":tlive_g4}, 
            "NR":{"E":evec_nr_nocap, "N":nhit_nr_nocap, "tlive":tlive_g4}}

############################################################################
#load simulated Capture data
#rcapture: expected capture rate, used to set livetime
# See calculations in: https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_n125:full_signal_fit#efficiencies_and_live_time
#load_frac: fraction of simulated events to load. Adjusts livetime appropriately
def load_simcap(file='data/v3_400k.pkl', rcapture=0.218, load_frac=1.0,verbose=False):
    if verbose:
        print('Loading (n,gamma) Data...')
    
    import pickle as pkl
    
    #load up some cascade simulated data
    with open(file,'rb') as readFile:
          cdata=pkl.load(readFile,encoding='latin1')

    #print(cdata.keys())
    if verbose:
        print(cdata['totalevents'])

    #Use only the events where the gamma escapes
    #Assume those where it doesn't escape end up at high Eee
    E_ng = cdata['E'][cdata['cEscape']]
    dE_ng = cdata['delE'][cdata['cEscape']]
    N = cdata['n'][cdata['cEscape']]
    
    nload=int(len(E_ng)*load_frac)
    
    E_ng = E_ng[:nload]
    dE_ng = dE_ng[:nload]
    N = N[:nload]
    
    #https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_n125:full_signal_fit&#efficiencies_and_live_time
    tlive_ng = load_frac*cdata['totalevents']/rcapture #[s]
    
    return {"E":E_ng, "dE":dE_ng, "N":N, "tlive":tlive_ng}

def load_simcap_old(lifetimes='fast', Ncascades='200k',verbose=False):
    if verbose:
        print('Loading (n,gamma) Data...')
    #'fast' or 'slow'
    #'2M' or '200k'
    
    import pickle as pkl
    
    #load up some cascade simulated data
    with open('data/normsi_{0}_{1}.pkl'.format(lifetimes, Ncascades),'rb') as readFile:
          cdata=pkl.load(readFile,encoding='latin1')

    #print(cdata.keys())
    if verbose:
        print(cdata['totalevents'])

    #Use only the events where the gamma escapes
    #Assume those where it doesn't escape end up at high Eee
    E_ng = cdata['E'][cdata['cEscape']]
    dE_ng = cdata['delE'][cdata['cEscape']]
    N = cdata['n'][cdata['cEscape']]
    
    #https://zzz.physics.umn.edu/cdms/doku.php?id=cdms:k100:run_summary:run_68:run_68_n125:full_signal_fit&#efficiencies_and_live_time
    tlive_ng = cdata['totalevents']/0.218 #[s]
    
    return {"E":E_ng, "dE":dE_ng, "N":N, "tlive":tlive_ng}