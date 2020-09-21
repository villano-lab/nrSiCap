import numpy as np
import dataPython as dp
import pandas as pd

#Plot previous yield measurements
#
#ax: The axis on which to plot
def plotOldYs(ax, **kwargs):
    #######################
    #Load data
    #######################

    #Damic
    ##########
    damic_data = dp.getXYdata_wXYerr('data/DAMIC_siyield_allerr.txt')

    #convert to numpy arrays
    damic_data['xx']= np.asarray(damic_data['xx'])*1000 #make units eV
    damic_data['yy']= np.asarray(damic_data['yy'])*1000 #make units eV
    damic_data['ex']= np.asarray(damic_data['ex'])*1000 #make units eV
    damic_data['ey']= np.asarray(damic_data['ey'])*1000 #make units eV

    #get the yield stuff
    damic_data['yy_yield'] = damic_data['yy']/damic_data['xx']
    damic_data['ey_yield'] = damic_data['yy_yield'] * np.sqrt((damic_data['ey']/damic_data['yy'])**2  + \
                                                              (damic_data['ex']/damic_data['xx'])**2)
    
    
    #Izraelevitch
    ##########
    #E_ion	dE_ion_sys	E_NR	dE_NR_stat	dE_NR_sys	Y	dY_stat	dY_sys
    izr_data = pd.read_csv("data/Izraelevitch_17.txt", sep='\t')
    #convert to numpy arrays, [eV]
    izr_data['E_NR']=np.asarray(izr_data['E_NR'])*1000
    izr_data['dE_NR_stat']=np.asarray(izr_data['dE_NR_stat'])*1000
    izr_data['dE_NR_sys']=np.asarray(izr_data['dE_NR_sys'])*1000  
    izr_data['Y']=np.asarray(izr_data['Y'])
    izr_data['dY_stat']=np.asarray(izr_data['dY_stat'])
    izr_data['dY_sys']=np.asarray(izr_data['dY_sys'])
    #Combine errors
    izr_data['dE_NR_tot']=np.sqrt(izr_data['dE_NR_stat']**2 + izr_data['dE_NR_sys']**2)
    izr_data['dY_tot']=np.sqrt(izr_data['dY_stat']**2 + izr_data['dY_sys']**2)
    
    #Dougherty
    ##########
    #Er	Y	dY-	dY+
    dough_data = pd.read_csv("data/Dougherty_91.txt", sep='\t')
    #convert to numpy arrays, [eV]
    dough_data['E_NR']=np.asarray(dough_data['Er'])*1000
    dough_data['Y']=np.asarray(dough_data['Y'])
    dough_data['dY-']=-1*np.asarray(dough_data['dY-'])#Stored as negative values
    dough_data['dY+']=np.asarray(dough_data['dY+'])

    
    #Gerbier
    ##########
    #Er	Y	dY-	dY+
    gerb_data = pd.read_csv("data/Gerbier_90.txt", sep='\t')
    #convert to numpy arrays, [eV]
    gerb_data['E_NR']=np.asarray(gerb_data['Er'])*1000
    gerb_data['Y']=np.asarray(gerb_data['Y'])
    gerb_data['dY-']=-1*np.asarray(gerb_data['dY-'])#Stored as negative values
    gerb_data['dY+']=np.asarray(gerb_data['dY+'])
    
    #Zecher
    ##########
    #Er	Y	dY-	dY+
    zech_data = pd.read_csv("data/Zecher_90.txt", sep='\t')
    #convert to numpy arrays, [eV]
    zech_data['E_NR']=np.asarray(zech_data['Er'])*1000
    zech_data['Y']=np.asarray(zech_data['Y'])
    zech_data['dY-']=-1*np.asarray(zech_data['dY-'])#Stored as negative values
    zech_data['dY+']=np.asarray(zech_data['dY+'])
    
    #Sattler
    ##########
    #Er	Y	dY-	dY+
    satt_data = pd.read_csv("data/Sattler_65.txt", sep='\t')
    #convert to numpy arrays, [eV]
    satt_data['E_NR']=np.asarray(satt_data['Er'])*1000
    satt_data['Y']=np.asarray(satt_data['Y'])
    satt_data['dY-']=-1*np.asarray(satt_data['dY-'])#Stored as negative values
    satt_data['dY+']=np.asarray(satt_data['dY+'])
    
    #CDMS II
    ##########
    #Er	dEr	Y	dY
    cdms2_data = pd.read_csv("data/CDMSII.txt", sep='\t')
    #convert to numpy arrays, [eV]
    cdms2_data['E_NR']=np.asarray(cdms2_data['Er'])*1000
    cdms2_data['dE_NR']=np.asarray(cdms2_data['dEr'])*1000
    cdms2_data['Y']=np.asarray(cdms2_data['Y'])
    cdms2_data['dY']=np.asarray(cdms2_data['dY'])

    
    #######################
    # Plot
    #######################
    ax.errorbar(damic_data['xx'], damic_data['yy_yield'], xerr=damic_data['ex'], yerr=damic_data['ey_yield'], \
                label='Damic', **kwargs)
    
    ax.errorbar(izr_data['E_NR'], izr_data['Y'], xerr=izr_data['dE_NR_tot'], yerr=izr_data['dY_tot'], \
                label='Izraelevitch', **kwargs)
    
    ax.errorbar(dough_data['E_NR'], dough_data['Y'], yerr=[dough_data['dY-'],dough_data['dY+']], \
                label='Dougherty', **kwargs)
    
    ax.errorbar(gerb_data['E_NR'], gerb_data['Y'], yerr=[gerb_data['dY-'],gerb_data['dY+']], \
                label='Gerbier', **kwargs)
        
    ax.errorbar(zech_data['E_NR'], zech_data['Y'], yerr=[zech_data['dY-'],zech_data['dY+']], \
                label='Zecher', **kwargs)
    
    ax.errorbar(satt_data['E_NR'], satt_data['Y'], yerr=[satt_data['dY-'],satt_data['dY+']], \
                label='Sattler', **kwargs)
    
    ax.errorbar(cdms2_data['E_NR'], cdms2_data['Y'], xerr=cdms2_data['dE_NR'], yerr=cdms2_data['dY'], \
                label='CDMS II', **kwargs)
    
    
#A function to load previous yield measurement data
#TODO: include more data...
def get_old_Y_data(label='izr'):
    if label.lower()=='izr':
        #Izraelevitch
        ##########
        #E_ion	dE_ion_sys	E_NR	dE_NR_stat	dE_NR_sys	Y	dY_stat	dY_sys
        izr_data = pd.read_csv("data/Izraelevitch_17.txt", sep='\t')
        #convert to numpy arrays, [eV]
        izr_data['E_NR']=np.asarray(izr_data['E_NR'])*1000
        izr_data['dE_NR_stat']=np.asarray(izr_data['dE_NR_stat'])*1000
        izr_data['dE_NR_sys']=np.asarray(izr_data['dE_NR_sys'])*1000  
        izr_data['Y']=np.asarray(izr_data['Y'])
        izr_data['dY_stat']=np.asarray(izr_data['dY_stat'])
        izr_data['dY_sys']=np.asarray(izr_data['dY_sys'])
        #Combine errors
        izr_data['dE_NR_tot']=np.sqrt(izr_data['dE_NR_stat']**2 + izr_data['dE_NR_sys']**2)
        izr_data['dY_tot']=np.sqrt(izr_data['dY_stat']**2 + izr_data['dY_sys']**2)
        
        return {'Enr':np.array(izr_data['E_NR']),'Y':np.array(izr_data['Y']),'dEnr':np.array(izr_data['dE_NR_tot']),'dY':np.array(izr_data['dY_tot'])}
    else:
        print("D'oh! We didn't include that one yet...")