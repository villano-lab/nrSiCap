#Likelihood and probabilty functions and things like that

import numpy as np
from numpy import log
from scipy import special
from scipy.special import factorial, gamma, loggamma

#Chisquared
def chisq(ydata,ypred,sd):
    return np.sum( ((ydata-ypred)/sd)**2 ) 

#Poisson likelihood of measuring k given expected mean of lambda
def pois_likelihood(k, lamb):
    return (lamb**k)*np.exp(-lamb)/gamma(k+1.)

#Poisson log-likelihood
#k: observed counts
#lamb: expected (model) counts
def ll_pois(k, lamb):   
    if np.sum(lamb<=0):
        return -np.inf
    
    return np.sum(k*log(lamb) - lamb - loggamma(k+1.))

#Normal log-likelihood, limit of Poisson for large lambda
#k: observed counts
#lamb: expected (model) counts
def ll_norm(k,lamb):
    if np.sum(lamb<=0):
        return -np.inf
    
    return np.sum(-0.5*log(2*np.pi*lamb) - (k-lamb)**2/(2*lamb))

#Log of flat prior functions
#theta: array of parameter values
#bounds: array of parameter bounds. shape should be len(theta)x2
def lp_flat(theta, bounds):
    #for itheta,ibounds in zip(theta,bounds):
    #    if not (ibounds[0] < itheta < ibounds[1]):
    #        return -np.inf
        
    #return 0.0
    
    if (np.array(bounds)[:,0]<=theta).all() and (theta<=np.array(bounds)[:,1]).all():
        #return 0.0
        return np.sum(-np.log(np.array(bounds)[:,1]-np.array(bounds)[:,0]))
    return -np.inf

#Log of normal prior distribution
#theta: parameter value(s)
#mu: parameter prior distribution mean(s)
#sigma: paramter prior distribution sigma(s)
def lp_norm(theta, mu, sigma):
    return np.sum(-0.5*((theta-mu)/sigma)**2 - log(sigma)-0.5*log(2*np.pi))


#Split-Normal log-likelihood
#k: test counts
#mu: mode counts
#sigma1: low-side sigma 
#sigma2: high-side sigma
def ll_SNorm(k,mu,sigma1,sigma2):
    if np.any(k<=0):
        return -np.inf
    #Initialize with ln(C)
    ll=0.5*log(2/np.pi)-log(sigma1+sigma2)
    #Handle the two cases
    ll[k<=mu]-=0.5*((k[k<=mu]-mu[k<=mu])/sigma1[k<=mu])**2
    ll[k>mu]-=0.5*((k[k>mu]-mu[k>mu])/sigma2[k>mu])**2
    #Return total log likelihood
    return np.sum(ll)

#Standard Normal probability density function
def fN(x):
    return (1/np.sqrt(2*np.pi))*np.exp(-(x**2)/2)
#Standard Normal CDF
def cN(k):
    return (1+special.erf(k/np.sqrt(2)))/2
#Quantile function of the standard normal distribution
def qN(q):
    return np.sqrt(2)*special.erfinv(2*q-1)

#Split-normal probability density function
def fSN(x,mu,sigma1,sigma2):
    x=np.asarray(x)
    prob=np.zeros_like(x)
    C=np.sqrt(2/np.pi)/(sigma1+sigma2)
    prob[x<=mu]=C*np.exp(-((x[x<=mu]-mu)/sigma1)**2/2)
    prob[x>mu]=C*np.exp(-((x[x>mu]-mu)/sigma2)**2/2)
    return prob

#CDF of Split-Normal (Eq 2.8 of Julio)
def cSN(k,mu,sigma1,sigma2):
    k=np.asarray(k)
    prob=np.zeros_like(k)
    C=np.sqrt(2/np.pi)/(sigma1+sigma2)
    prob[k<=mu]=C*np.sqrt(2*np.pi)*sigma1*cN((k[k<=mu]-mu)/sigma1)
    prob[k>mu]=1-C*np.sqrt(2*np.pi)*sigma2*(1-cN((k[k>mu]-mu)/sigma2))
    return prob

#Quantile function of the Split-Normal distribution (Eq 2.9 of Julio)
def qSN(q,mu,sigma1,sigma2):
    q=np.asarray(q)
    C=np.sqrt(2/np.pi)/(sigma1+sigma2)
    p=cSN(mu,mu,sigma1,sigma2)
    k=np.zeros_like(q)
    k[q<=p]=mu+sigma1*qN(q[q<=p]/(C*sigma1*np.sqrt(2*np.pi)))
    k[q>p]=mu+sigma2*qN((q[q>p]+C*sigma2*np.sqrt(2*np.pi)-1)/(C*sigma2*np.sqrt(2*np.pi)))
    return k

#For a given set of median and upper and lower quantile points, find the best fit split gaussian to them
#xs: points
#qs: quantiles, e.g. [0.15865,0.5,0.84135] for median and +/- 1-sigma
#w: array of weighting factors for target points. Default = uniform
#mode: Can perscribe a mode exactly
#bounds: bounds of mu,sigma1,sigma2 values in the form ((mu_low,muhi),(sigma1_low,sigma1_hi),(sigma2_low,sigma2_hi))
#   Values can be None. Recommend using small lower values to keep them non-zero.
#Returns: best fit mu,sigma1,sigma2
from scipy.optimize import minimize
def getSNpars(xs,qs,w=None,mode=None,bounds=((1e-12,None),(1e-12,None),(1e-12,None))):
    xs=np.asarray(xs)
    qs=np.asarray(qs)
    if w is None:
        w=np.ones_like(xs)
    else:
        w=np.asarray(w)
        
    if mode is None:
        #theta = [mode,sigma1,sigma2]
        chisq = lambda theta: np.sum(w*(xs - qSN(qs,*theta))**2)
        #Crude estimate of mode and sigmas to initialize fit
        pfit = np.poly1d(np.polyfit(qs,xs,deg=1))
        #plot(pfit(qs),qs)
        x0=pfit([0.5,0.15865,0.84135])
        theta0=np.array([x0[0],x0[0]-x0[1],x0[2]-x0[0]])

        res=minimize(chisq,theta0,bounds=bounds)
        return res['x']
        
    else:
        #theta = [sigma1,sigma2], mode is fixed
        chisq = lambda theta: np.sum(w*( xs - qSN(qs,mode,*theta))**2)
        #We know the mode exactly, just need sigmas
        pfit = np.poly1d(np.polyfit([0.5,*qs],[mode,*xs],deg=1))
        x0=pfit([0.15865,0.84135])
        theta0=np.array([mode-x0[0],x0[1]-mode])
        
        res=minimize(chisq,theta0,bounds=bounds[1:])
        return np.insert(res['x'],0,mode)
    
#Fit Split-Normals to an array of measured and sigma values
def getSNparsArray(modes,sigma_up, sigma_down):
    popt=[]
    for i in np.arange(len(modes)):
        popt.append(getSNpars([modes[i]-sigma_down[i],modes[i]+sigma_up[i]],[0.15865,0.84135],w=[1,1],mode=modes[i]))
    return np.array(popt)