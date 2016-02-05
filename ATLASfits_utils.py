#--ATLASfits_utils.py - Version 1 - 04/02/2016
#--Author: Bradley J Kavanagh
#--Summary: Auxiliary functions to be used in ATLASfits.py
#--and the associated plotting routines
#--Note: Requires emcee (http://dan.iel.fm/emcee/current/)
#--Please report any problems to: bradley.kavanagh@lpthe.jussieu.fr

import sys
import numpy as np
from scipy.special import erf
from scipy.integrate import quad
from scipy.stats import chi2
import emcee


#---Options--
#Set this to 1 if you want to exclude empty bins from the analysis
exclude_bins = 0

#---Initialisation (load ATLAS data etc.)----

sqrts = 13000

#Load ATLAS data
data = np.loadtxt("ATLASdata1.txt")[:,1]
nbins = data.size

#Calculate bin edges
m_edges = np.linspace(150,150+40*(nbins),nbins+1)
m_centres = m_edges[:-1] + 20

n_emptybins = np.sum(data < 1)

used_bins = nbins
if (exclude_bins):
    used_bins = used_bins - n_emptybins

if (exclude_bins):
    print "   Using", used_bins, " bins [excluding empty bins]..."
else:
    print "   Using", used_bins, " bins [including empty bins]..."

#Define number of parameters in various signal and BG models
Nbgparams = [2, 3, 4, 3, 4, 5]
Nsigparams = [0,2,3]


#---Signal functions-------------

#Signal distribution for a wide resonance
def f_wide(x, M, alpha):
    gamma = M**2*np.sqrt(1 + (alpha)**2)
    k = 2*np.sqrt(2)*(M**2)*(alpha)*gamma/(np.pi*np.sqrt(M**2 + gamma))
    return k*1.0/((x**2 - M**2)**2 + M**4*alpha**2)

#Signal distribution for a narrow resonance
def f_NWA(x, M, sigma):
    return np.exp(-0.5*((x - M)**2)/sigma**2)/(np.sqrt(2*np.pi)*sigma)

#Diphoton invariant mass resolution
def resolution(M):
    return 6.11e-3*M + 0.7778
    

#----Background functions---------

#k = 0
def f_ATLAS0(x, b, a0):
    return ((1-(x/sqrts)**(1.0/3.0))**b)*((x/sqrts)**a0)
    
#k = 1
def f_ATLAS1(x, b, a0, a1):
    return ((1-(x/sqrts)**(1.0/3.0))**b)*((x/sqrts)**(a0  + a1*np.log(x/sqrts)))

#k = 2
def f_ATLAS2(x, b, a0, a1, a2):
    return ((1-(x/sqrts)**(1.0/3.0))**b)*((x/sqrts)**(a0  + a1*np.log(x/sqrts) + a2*np.log(x/sqrts)**2))
    

#-------Priors--------

def lnprior_BG(x, BG_type):
    for i in range(Nbgparams[BG_type]):
        if (np.abs(x[i]) > 25):
            return -np.inf
    return 0

def lnprior_sig(logA, M , alpha = 5.0):
    if ((logA < -2) or (logA > 2)):
        return -np.inf
    if ((M < 700)or(M > 800)):
        return -np.inf
    if ((alpha < 1.0)or(alpha > 10.0)):
        return -np.inf
    return 0


#----Calculate log-likelihood--------

def lnprob(x, BG_type, signal_type):
    
    #Check against BG priors
    lp = lnprior_BG(x, BG_type)
    if not np.isfinite(lp):
            return -np.inf
            
    #Specify number of parameters
    Nbg = Nbgparams[BG_type]
    Nsig = Nsigparams[signal_type]
    ndim = Nbg + Nsig
    

    ind = 0
    BG_norm = 1.0
    #Unpack background variables
    if (BG_type in [3, 4, 5]): #Free normalisation
        BG_norm = 10**x[0]
        ind = 1
    if (BG_type in [0, 3]):
        b = x[0+ind]
        a0 = x[1+ind]
    elif (BG_type in [1, 4]):
        b = x[0+ind]
        a0 = x[1+ind]
        a1 = x[2+ind]
    elif (BG_type in [2, 5]):
        b = x[0+ind]
        a0 = x[1+ind]
        a1 = x[2+ind]
        a2 = x[3+ind]
    
    #Unpack and calculate signal variables
    # M - resonance mass
    # logA - log_10(Number of signal events)
    # alpha - width/mass
    if (signal_type > 0):
        M = x[Nbg]
        if (signal_type == 1): #Narrow resonance
            logA = x[Nbg+1]
            lp_sig = lnprior_sig(logA, M)
        elif (signal_type == 2): #Wide resonance
            alpha = x[Nbg+1]
            logA = x[Nbg+2] 
            lp_sig = lnprior_sig(logA, M, alpha)
            
        if not np.isfinite(lp_sig):
            return -np.inf
    
    #Calculate the likelihood
    lnL = 0
    
    #Sum over bins
    for i in range(nbins):
        No = data[i]
        Ne = 0

        #Exclude bins if required
        if ((exclude_bins)and(No < 0.5)):
            lnL = lnL + 0
        else:
            #Calculate background contribution
            if (BG_type in [0, 3]): #k = 0
                Ne = BG_norm*quad(f_ATLAS0, m_edges[i], m_edges[i+1], args=(b, a0))[0]/40.0
            elif (BG_type in [1, 4]): #k = 1
                Ne = BG_norm*quad(f_ATLAS1, m_edges[i], m_edges[i+1], args=(b, a0, a1))[0]/40.0
            elif (BG_type in [2, 5]): #k = 2
                Ne = BG_norm*quad(f_ATLAS2, m_edges[i], m_edges[i+1], args=(b, a0, a1, a2))[0]/40.0
                
            #Calculate signal contribution 
            N_S = 0
            if (signal_type == 1): #Narrow resonance
                N_S = (10**logA)*quad(f_NWA, m_edges[i], m_edges[i+1], args=(M,resolution(M)))[0]

            elif (signal_type == 2): #Wide resonance
                N_S = (10**logA)*quad(f_wide, m_edges[i], m_edges[i+1], args=(M,alpha/100.0))[0]

            Ne = Ne+N_S
            
            #Check that everything is fine
            if (Ne < 1e-10):
                Ne = 1e-10      
            if (Ne < 0):
                print "ATLASfits.py (near line 200): You messed up..."
            if not np.isfinite(Ne):
                return -np.inf
        
            #Poisson likelihood
            lnL = lnL - Ne + No*np.log(Ne)

    return lnL


#----Calculate best fit likelihood and parameter values----

def getBestFit(BG_type, signal_type):
    #emcee parameters
    nwalkers, nsamps = 20, 1000
    
    #Specify number of parameters
    Nbg = Nbgparams[BG_type]
    Nsig = Nsigparams[signal_type]
    ndim = Nbg + Nsig
    
    #Calculate the constant term in the likelihood
    lnNo = 0
    for i in range(nbins):
        No = data[i]
        if (No > 0):
            for j in range(int(No)):
                lnNo += np.log(j+1)
      
    #Specify the likelihood function  
    lnp = lambda x: lnprob(x, BG_type, signal_type)
    
    #Specify initial parameter values x0
    #These are the best-fit BG-only parameters I've found
    if (BG_type == 0): #Without free normalisation
        x0 = [12.77024061 , -2.68784895]
    elif (BG_type == 1):
        x0 = [ 9.58916916 ,-1.72585052,  0.17706063]
    elif (BG_type == 2):
        x0 =  [2.4173727 ,  3.57063531 , 2.19234225 , 0.20673052]
    
    elif (BG_type == 3): #With free normalisation
        x0 = [-3.07349579 , 4.19769785, -3.79149089]
    elif (BG_type == 4):
        x0 = [-15.81710377 ,-17.60126389, -12.51186814 , -0.7620349]
    elif (BG_type == 5):
        x0 = [13.89459469 , 20.45122052  ,17.37188926  , 4.66081521 ,  0.37492859]

    #Append signal parameters
    if (signal_type == 1): #NWA
        x0_s = [750.0, 1.0]
    elif (signal_type == 2): #Wide
        x0_s = [750.0, 5.0, 1.0]
    
    if (signal_type > 0):
        x0 = np.append(x0, x0_s)
    
    p0 = [(x0 + np.random.rand(ndim)*1e0) for i in range(nwalkers)]

    #Run the MCMC sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnp, args=())
    sampler.run_mcmc(p0, nsamps)

    #Extract samples and likelihoods
    samples = sampler.flatchain 
    likes = sampler.flatlnprobability - lnNo

    #Return maximum log-likelihood (lnL) and best fit points
    lnL = likes[np.argmax(likes)]
    bf = samples[np.argmax(likes)]
    return lnL, bf
    
#-----------------------------------