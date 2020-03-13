#NR ionization yield class
import numpy as np
from scipy.special import erf
import sys
sys.path.append('../python/')
import damic_y as dy
#see the link below for the reason for the filtrations
#https://stackoverflow.com/questions/40845304/runtimewarning-numpy-dtype-size-changed-may-indicate-binary-incompatibility
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")

class Yield:
    def __init__(self, model, pars):
        self.models={'Lind': 'Lindhard','Chav': 'Chavarria', 'Sor': 'Sorenson', 'Damic': 'Extrapolated Damic model', 'AC': 'Adiabatic Correction'}
        self.model_npar={'Lind': 1,'Chav': 2, 'Sor': 2, 'Damic': 0, 'AC': 2}
        self.set_model(model)
        self.set_pars(pars)
        
    def set_model(self,model):
        if model in self.models.keys():
            self.model=model
        else:
            print('Error: '+str(model)+' is not a valid Yield model. Please select from: '+str(self.models.keys()))
    
    def set_pars(self,pars):
        if len(pars) == self.model_npar[self.model]:
            self.pars=pars
            self.npars=len(self.pars)
        else:
            print('Error: '+str(self.model)+' yield model takes '+str(self.model_npar[self.model])+' parameter(s), but '+str(len(pars))+' are given.')
                 
    def calc(self,Er): #expects NR energy in [eV]
        if len(self.pars) != self.model_npar[self.model]:
            print('Error: '+str(self.model)+' yield model takes '+str(self.model_npar[self.model])+' parameter(s), but '+str(len(self.pars))+' are given.')
            return None
            
        if self.model=='Lind':
                return yLind(Er, self.pars[0])        
        elif self.model=='Chav':
                return yChav(Er, self.pars[0], self.pars[1]) 
        elif self.model=='Sor':
                return ySor(Er, self.pars[0], self.pars[1])
        elif self.model=='Damic':
                return yDamic(Er)     
        elif self.model=='AC':
                return yAC(Er,self.pars[0], self.pars[1])
        else:
            print('Error: '+str(self.model)+' yield model not defined')
            return None
        

#Lindhard yield
#http://gymarkiv.sdu.dk/MFM/kdvs/mfm%2030-39/mfm-33-10.pdf
#k = 0.133Z^(2/3)A^(−1/2)
#This is only decent for eps>0.01 => Er>400 eV

# Er: Nuclear recoil energy [eV]
# k: Lindhard k factor [unitless]
def yLind(Er, k):
    Z=14.
    eps = 11.5*Er/1000*Z**(-7./3)
    g = 3.*eps**0.15 + 0.7*eps**0.6 + eps
    return (k*g)/(1+k*g)

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
    Z=14.
    eps = 11.5*Er/1000*Z**(-7./3)
    
    y=yLind(Er,k)
    
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


#a spline extrapolation to DAMIC data
yDamic = np.vectorize(dy.getDAMICy())