import numpy as np
import dataPython as dp
import pandas as pd
import os

#Reference to this file
def package_path(*paths, package_directory=os.path.dirname(os.path.abspath(__file__))):
    return os.path.join(package_directory, *paths)

#Plot previous yield measurements
#
#ax: The axis on which to plot
def plotOldYs(ax, datasets=['chav','izr','dough','gerb','zech','satt','agnese'],
              labels=['Chavarria','Izraelevitch','Dougherty','Gerbier','Zecher','Sattler','Agnese'], **kwargs):

    for ds,label in zip(datasets,labels):
        #Load data
        data=get_old_Y_data(label=ds)
        
        Enr=data['Enr']
        Y=data['Y']
        
        if 'dEnr' in data.keys():
            dEnr=data['dEnr']
        else:
            dEnr=None
            
        if 'dY' in data.keys():
            dY=data['dY']
        elif 'dY-' in data.keys():
            dY=np.array([data['dY-'],data['dY+']])
        else:
            dY=None
            
        #Plot
        ax.errorbar(Enr,Y, xerr=dEnr, yerr=dY, label=label, **kwargs)


def plotOldYs_noSat(ax, datasets=['chav','izr','dough','gerb','zech','agnese'],
              labels=['Chavarria','Izraelevitch','Dougherty','Gerbier','Zecher','Agnese'], **kwargs):

    for ds,label in zip(datasets,labels):
        #Load data
        data=get_old_Y_data(label=ds)
        
        Enr=data['Enr']
        Y=data['Y']
        
        if 'dEnr' in data.keys():
            dEnr=data['dEnr']
        else:
            dEnr=None
            
        if 'dY' in data.keys():
            dY=data['dY']
        elif 'dY-' in data.keys():
            dY=np.array([data['dY-'],data['dY+']])
        else:
            dY=None
            
        #Plot
        ax.errorbar(Enr,Y, xerr=dEnr, yerr=dY, label=label, **kwargs)
        
        
#A function to load previous yield measurement data
def get_old_Y_data(label='izr'):
    
    if 'chav' in label.lower():
        #Damic
        #https://arxiv.org/pdf/1608.00957.pdf
        ##########
        #Ee	Er	dEr_stat	dEr_sys
        chav_data = pd.read_csv(package_path()+"/../data/Chavarria_16.txt", sep='\t')
        #convert to numpy arrays, [eV]
        chav_data['Ee']=np.asarray(chav_data['Ee'])*1000
        chav_data['Er']=np.asarray(chav_data['Er'])*1000
        chav_data['dEr_stat']=np.asarray(chav_data['dEr_stat'])*1000
        chav_data['dEr_sys']=np.asarray(chav_data['dEr_sys'])*1000
        
        chav_data['dEr_tot']=np.sqrt(chav_data['dEr_stat']**2 + chav_data['dEr_sys']**2)
        chav_data['Y']=chav_data['Ee']/chav_data['Er']
        chav_data['dY']=chav_data['Y']*(chav_data['dEr_tot']/chav_data['Er'])

        
        return {'Enr':np.array(chav_data['Er']),
                'Y':np.array(chav_data['Y']),
                'dEnr':np.array(chav_data['dEr_tot']),
                'dY':np.array(chav_data['dY'])}

    
    elif 'izr' in label.lower():
        #Izraelevitch
        #https://arxiv.org/pdf/1702.00873.pdf
        ##########
        #E_ion	dE_ion_sys	E_NR	dE_NR_stat	dE_NR_sys	Y	dY_stat	dY_sys
        izr_data = pd.read_csv(package_path()+"/../data/Izraelevitch_17.txt", sep='\t')
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
        
        return {'Enr':np.array(izr_data['E_NR']),
                'Y':np.array(izr_data['Y']),
                'dEnr':np.array(izr_data['dE_NR_tot']),
                'dY':np.array(izr_data['dY_tot'])}

    elif 'dough' in label.lower():
        #Dougherty
        #10.1103/PhysRevA.45.2104
        ##########
        #Er	dEr	Y	dY
        dough_data = pd.read_csv(package_path()+"/../data/Dougherty_91.txt", sep='\t')
        #convert to numpy arrays, [eV]
        dough_data['Er']=np.asarray(dough_data['Er'])*1000
        dough_data['dEr']=np.asarray(dough_data['dEr'])*1000
        dough_data['Y']=np.asarray(dough_data['Y'])
        dough_data['dY']=np.asarray(dough_data['dY'])

        return {'Enr':np.array(dough_data['Er']),
                'Y':np.array(dough_data['Y']),
                'dEnr':np.array(dough_data['dEr']),
                'dY':np.array(dough_data['dY'])}
    
    elif 'gerb' in label.lower():
        #Gerbier
        #10.1103/PhysRevD.42.3211
        ##########
        #Er	dEr	Y	dY
        gerb_data = pd.read_csv(package_path()+"/../data/Gerbier_90.txt", sep='\t')
        #convert to numpy arrays, [eV]
        gerb_data['Er']=np.asarray(gerb_data['Er'])*1000
        gerb_data['dEr']=np.asarray(gerb_data['dEr'])*1000
        gerb_data['Y']=np.asarray(gerb_data['Y'])
        gerb_data['dY']=np.asarray(gerb_data['dY'])

        return {'Enr':np.array(gerb_data['Er']),
                'Y':np.array(gerb_data['Y']),
                'dEnr':np.array(gerb_data['dEr']),
                'dY':np.array(gerb_data['dY'])}
    
    elif 'zech' in label.lower():
        #Zecher
        #10.1103/PhysRevA.41.4058
        ##########
        #Er	dEr	Y	dY-	dY+
        zech_data = pd.read_csv(package_path()+"/../data/Zecher_90.txt", sep='\t')
        #convert to numpy arrays, [eV]
        zech_data['Er']=np.asarray(zech_data['Er'])*1000
        zech_data['dEr']=np.asarray(zech_data['dEr'])*1000
        zech_data['Y']=np.asarray(zech_data['Y'])
        zech_data['dY-']=np.asarray(zech_data['dY-'])#Stored as negative values
        zech_data['dY+']=np.asarray(zech_data['dY+'])
        
        return {'Enr':np.array(zech_data['Er']),
                'Y':np.array(zech_data['Y']),
                'dEnr':np.array(zech_data['dEr']),
                'dY-':np.array(zech_data['dY-']),
                'dY+':np.array(zech_data['dY+'])}
    
    elif 'satt' in label.lower():
        #Sattler
        #10.1103/PhysRev.138.A1815
        ##########
        #N.B. I don't know if these are right. I think there should be a lot more values and
        #  there should be a decently large Er uncertainty that's not in their table I
        #Er	Y	dY-	dY+
        satt_data = pd.read_csv(package_path()+"/../data/Sattler_65.txt", sep='\t')
        #convert to numpy arrays, [eV]
        satt_data['Er']=np.asarray(satt_data['Er'])*1000
        satt_data['Y']=np.asarray(satt_data['Y'])
        satt_data['dY-']=np.asarray(satt_data['dY-'])
        satt_data['dY+']=np.asarray(satt_data['dY+'])
        
        return {'Enr':np.array(satt_data['Er']),
                'Y':np.array(satt_data['Y']),
                'dY-':np.array(satt_data['dY-']),
                'dY+':np.array(satt_data['dY+'])}
    
    elif 'agnese' in label.lower():
        #CDMS-II
        #10.1016/j.nima.2018.07.028
        #Er	dEr	Y	dY
        cdms2_data = pd.read_csv(package_path()+"/../data/CDMSII.txt", sep='\t')
        #convert to numpy arrays, [eV]
        cdms2_data['Er']=np.asarray(cdms2_data['Er'])*1000
        cdms2_data['dEr']=np.asarray(cdms2_data['dEr'])*1000
        cdms2_data['Y']=np.asarray(cdms2_data['Y'])
        cdms2_data['dY']=np.asarray(cdms2_data['dY'])
        
        return {'Enr':np.array(cdms2_data['Er']),
                'Y':np.array(cdms2_data['Y']),
                'dEnr':np.array(cdms2_data['dEr']),
                'dY':np.array(cdms2_data['dY'])}

    else:
        print("D'oh! We didn't include that one yet...")
        return None