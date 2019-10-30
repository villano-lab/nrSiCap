#NR ionization yield class
import numpy as np
from scipy.special import erf
import sys
sys.path.append('../python/')
import damic_y as dy 

class Yield:
    def __init__(self, model, pars):
        self.model = model
        self.pars=pars
        self.models=['Lind: Lindhard','Chav: Chavarria', 'Sor: Sorenson', 'Damic: Extrapolated Damic model']

    def calc(self,Er): #expects NR energy in [eV]
        
        if self.model=='Lind':
            if len(self.pars)!=1:
                print('Error: '+str(self.model)+' yield model takes 1 parameter, but '+len(self.pars)+' are given.')
                return None
            else:
                return yLind(Er, self.pars[0])
            
        elif self.model=='Chav':
            if len(self.pars)!=2:
                print('Error: '+str(self.model)+' yield model takes 2 parameters, but '+len(self.pars)+' are given.')
                return None
            else:
                return yChav(Er, self.pars[0], self.pars[1])
            
        elif self.model=='Sor':
            if len(self.pars)!=2:
                print('Error: '+str(self.model)+' yield model takes 2 parameters, but '+len(self.pars)+' are given.')
                return None
            else:
                return ySor(Er, self.pars[0], self.pars[1])
        
        elif self.model=='Damic':
            if len(self.pars)!=0:
                print('Error: '+str(self.model)+' yield model takes 0 parameterd, but '+len(self.pars)+' are given.')
                return None
            else:
                return yDamic(Er)
        else:
            print('Error: '+str(self.model)+' yield model not defined')
            return None
        

#Lindhard yield
#http://gymarkiv.sdu.dk/MFM/kdvs/mfm%2030-39/mfm-33-10.pdf
#k = 0.133Z^(2/3)A^(−1/2)
#This is only decent for eps>0.01 => Er>400 eV
def yLind(Er, k):
    Z=14.
    eps = 11.5*Er/1000*Z**(-7./3)
    g = 3.*eps**0.15 + 0.7*eps**0.6 + eps
    return (k*g)/(1+k*g)

#Lindhard w/ Chavarria tweak
#https://arxiv.org/pdf/1803.02903.pdf
def yChav(Er,k,a):
    return 1/(1/(a*Er)+1/yLind(Er,k))

#Sorenson: Lindhard + constant
#https://journals.aps.org/prd/pdf/10.1103/PhysRevD.91.083509
def ySor(Er,k,q):
    Z=14.
    eps = 11.5*Er/1000*Z**(-7./3)
    #This can go negative for large q/eps
    return np.maximum(0,yLind(Er,k)-q/eps)

#a spline extrapolation to DAMIC data
yDamic = np.vectorize(dy.getDAMICy())