#NR ionization yield class
#see the link below for the reason for the filtrations
#https://stackoverflow.com/questions/40845304/runtimewarning-numpy-dtype-size-changed-may-indicate-binary-incompatibility
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

import numpy as np
from scipy.special import erf
from scipy.interpolate import CubicSpline, PchipInterpolator
import sys
sys.path.append('../python/')
import damic_y as dy


#Conversion from eV to unitless energy, epsilon as defined in Lindhard yield
eVTOeps = 11.5/1000*14**(-7./3)

################################################################################
class Yield:
    def __init__(self, model, pars):
        self.models={'Lind': 'Lindhard','Chav': 'Chavarria', 'Sor': 'Sorenson', 'Damic': 'Extrapolated Damic model', 'AC': 'Adiabatic Correction', 'pchip':'Lindhard+PCHIP', 'Shexp':'Lindhard+shelf+exp', 'Pol3':'3-degree Polynomial', 'Pol4':'4-degree Polynomial'}
        self.model_npar={'Lind': 1,'Chav': 2, 'Sor': 2, 'Damic': 0, 'AC': 2, 'pchip':5, 'Shexp':5, 'Pol3':3, 'Pol4':4}
        self.set_model(model)
        self.set_pars(pars)
        
########################################      
    def set_model(self,model):
        if model in self.models.keys():
            self.model=model
        else:
            print('Error: '+str(model)+' is not a valid Yield model. Please select from: '+str(self.models.keys()))
            
########################################    
    def set_pars(self,pars):
        if len(pars) == self.model_npar[self.model]:
            self.pars=pars
            self.npars=len(self.pars)
        else:
            print('Error: '+str(self.model)+' yield model takes '+str(self.model_npar[self.model])+' parameter(s), but '+str(len(pars))+' are given.')

########################################
    #Precalculate anything that we can
    def solve(self):
        if self.model=='pchip':
            (Er0,Er1,Er2,f1)=self.pars[1:]
            Er_ad=np.array([Er0,Er1,Er2])
            f_ad=np.array([0,f1,1])
            #Check if ascending
            if not np.all(np.diff(Er_ad)>=0):
                #print('Error: pchip yield model requires ascending energy values, but received ',Er_ad)
                return False
            
            if not np.all(np.diff(f_ad)>=0):
                #print('Error: pchip yield model requires ascending fraction values, but received ',f_ad)
                return False
            
            self.f_pchip = PchipInterpolator(Er_ad, f_ad)
            return True
        else:
            return True
            
########################################                
    def calc(self,Er): #expects NR energy in [eV]
        if len(self.pars) != self.model_npar[self.model]:
            print('Error: '+str(self.model)+' yield model takes '+str(self.model_npar[self.model])+' parameter(s), but '+str(len(self.pars))+' are given.')
            return None
            
        if self.model=='Lind':
            return yLind(Er,*(self.pars))        
        elif self.model=='Chav':
            return yChav(Er, *(self.pars))
        elif self.model=='Sor':
            return ySor(Er, *(self.pars))
        elif self.model=='Damic':
            return yDamic(Er)     
        elif self.model=='AC':
            return yAC(Er,*(self.pars))
        elif self.model=='pchip':
            return yLind_pchip(Er,*(self.pars),self.f_pchip)
        elif self.model=='Shexp':
            return yLind_shelf_exp(Er,*(self.pars))
        elif self.model=='Pol3':
            return yPol3(Er,*(self.pars))
        elif self.model=='Pol4':
            return yPol4(Er,*(self.pars))
        else:
            print('Error: '+str(self.model)+' yield model not defined')
            return None
        
################################################################################

#Lindhard yield
#http://gymarkiv.sdu.dk/MFM/kdvs/mfm%2030-39/mfm-33-10.pdf
#k = 0.133Z^(2/3)A^(−1/2)
#This is only decent for eps>0.01 => Er>400 eV

# Er: Nuclear recoil energy [eV]
# k: Lindhard k factor [unitless]
def yLind(Er, k):
    kg = k*(3.*(Er*eVTOeps)**0.15 + 0.7*(Er*eVTOeps)**0.6 + (Er*eVTOeps))
    return (kg)/(1+kg)

#Lindhard w/ Chavarria tweak
#https://arxiv.org/pdf/1803.02903.pdf

# Er: Nuclear recoil energy [eV]
# k: Lindhard k factor [unitless]
# ainv: Inverse Chavaria cutoff factor [eV]
def yChav(Er,k,ainv):
    #return 1/(1/(a*Er/1000)+1/yLind(Er,k))
    
    #Avoid divide by 0's
    #A=a*Er/1000
    #B=yLind(Er,k)
    #y=np.zeros_like(A*B)
    #nz=((A!=0) & (B!=0))
    #y[nz] = A[nz]*B[nz]/(A[nz]+B[nz])
    #return y

    #Switch to ainv = 1/a parameterization
    yL=yLind(Er,k)
    y=Er*yL #Numerator
    y_denom=Er+yL*ainv
    
    #Handle lists and floats differently
    if isinstance(y_denom,(list, tuple, np.ndarray)):
        zero=(y_denom==0)
        y[~zero]/=y_denom[~zero]
        y[zero]=0
    else:
        if y_denom==0:
            y=0
        else:
            y/=y_denom
            
    return y
    
#Sorenson: Lindhard + constant
#https://journals.aps.org/prd/pdf/10.1103/PhysRevD.91.083509
#Er: Recoil energy [eV]
#k: Lindhard k parameter
#q: Sorenson cutoff paramter in unitless energy [eps]
def ySor(Er,k,q):
    y=yLind(Er,k)
    
    eps = Er*eVTOeps
    
    #Handle lists and floats differently
    if isinstance(eps,(list, tuple, np.ndarray)):
        #Avoid divide by 0
        y[eps>0]-=q/eps[eps>0]
        #Zero non-physical values
        y[y<0]=0
    else:
        if eps>0:
            y-=q/eps
        if y<0:
            y=0
    
    return y

#Lindhard with Adiabatic Correction
#https://journals.aps.org/prd/pdf/10.1103/PhysRevD.94.122003
#Er: Recoil Energy [eV]
#k: Lindhard k parameter
#xi: adiabatic energy scale factor
def yAC(Er,k,xi):
    return (1-np.exp(-Er/xi))*yLind(Er,k)

#Lindhard with monotonic cubic spline curve at low energy
#Spline takes over from Lind below energy Er2
#Yield goes to 0 at Er0
#Intermediate point at (Er1,f1) sets Y(Er1)=f1*Y_lind(Er2)
# i.e. how much Yield decreased since spline took over.
#
#This one requires pre-caculating the spline
#Parameters you need to set are k,Er0,Er1,Er2,f1
def yLind_pchip(Er,k,Er0,Er1,Er2,f1,f_pchip):
    y=yLind(Er, k)
    ind_E0_E2=(Er>=Er0) & (Er<Er2)
    
    y[ind_E0_E2]=f_pchip(Er[ind_E0_E2])*yLind(Er2, k)
    y[Er<Er0]=0
    
    return y

#Homebrew Lindhard with shelf and exp cutoff at low energy
#k: Lindhard k
#Yshelf: Value of flat yield shelf
#Ec: cutoff energy at which Y->0 [eV]
#dE: width of cutoff transition [eV]
#alpha: smoothing parameter between shelf and lindhard
def yLind_shelf_exp(Er,k,Yshelf,Ec,dE,alpha):
    yL=yLind(Er, k)

    return (Er>Ec)*(1-np.exp(-(Er-Ec)/dE))*np.power(np.power(yL,alpha)+np.power(Yshelf,alpha),1/alpha)

#a spline extrapolation to DAMIC data
yDamic = np.vectorize(dy.getDAMICy())

#Simple polynomials
def yPol3(E,p0,p1,p2):
    poly = np.poly1d((p2,p1,p0))(E)
    poly[poly<0]=0
    return poly

def yPol4(E,p0,p1,p2,p3):
    poly = np.poly1d((p3,p2,p1,p0))(E)
    poly[poly<0]=0
    return poly