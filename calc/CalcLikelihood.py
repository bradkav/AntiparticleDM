import sys
from scipy.stats import chi2, norm
#import matplotlib.pyplot as pl
import numpy as np

#Sampling paramaters
logc_min = -10
logc_max = -6

#Initial half-width of the grid
d0 = (logc_max - logc_min)/2.0

#Number of refinement steps
Ndvals = 10

#Calculate half-width of grid at each step
#   -->Shrink with each step
r = 1.5
dvals = d0/np.append((1 + np.exp(r*np.linspace(-4, 2,Ndvals-1))), 50.0)

#Should we plot likelihoods at intermediate steps?
PLOTLIKE = False

#----Functions----

def CalcLike_grid(mx, expts, Ngrid = 100, maj = False, refine=False):
    
    N_expt = len(expts)

    #Initial grid of parameters
    cp_list = np.logspace(logc_min, logc_max, Ngrid)
    cn_list = np.logspace(logc_min, logc_max, Ngrid)
    if (maj):
        f_list = np.asarray([-1.0, 1.0])
    else:
        f_list = np.linspace(-1.0,1.0, Ngrid)
    Nfvals = len(f_list)

    #Create 3-D mesh
    (CP, CN, F) = np.meshgrid(cp_list, cn_list, f_list, indexing='ij')

    #Calculate likelihood summing over experiments
    full_like = 0.0
    for expt in expts:
        No = len(expt.events)
        like = 0.0
        A = np.zeros((expt.N_iso, Ngrid, Ngrid, Nfvals))
    
        #Calculate a rescaling function which depends on the couplings
        for i in range(expt.N_iso):
            A[i, :, :, :] = 2.0*((CP*expt.N_p[i] + CN*expt.N_n[i])**2\
                     + 2.0*CP*CN*(F-1)*expt.N_p[i]*expt.N_n[i])

        A = np.clip(A, 1e-50, 1e50)

        #For a single element, use Eq. (A1) for likelihood
        if (expt.N_iso == 1): 
            like = -A[0,:,:,:]*expt.Ne_list
            like += expt.eventlike + No*np.log(A[0,:,:,:])
            
        #Otherwise we use Eq. (A2)
        else:   
            Ne_tot = np.einsum('ijkl,i', A, expt.Ne_list)
            like = -Ne_tot + No*np.log(Ne_tot)
            Nek = np.einsum('ijkl,i->ijkl', A, expt.Ne_list)/Ne_tot
            like += expt.eventlike_interp.ev(Nek[0,:,:,:],Nek[1, :,:,:])

        full_like += like

    if (maj):
        #Get best fit for Majorana-like case
        ind_maj_minus = np.argmax(full_like[:,:,0].flatten())
        cpmax_maj_minus = CP[:,:,0].flatten()[ind_maj_minus]
        cnmax_maj_minus = CN[:,:,0].flatten()[ind_maj_minus]
    
        ind_maj_plus = np.argmax(full_like[:,:,-1].flatten())
        cpmax_maj_plus = CP[:,:,-1].flatten()[ind_maj_plus]
        cnmax_maj_plus = CN[:,:,-1].flatten()[ind_maj_plus]
   
    else:
        #Get best fit for Dirac-like case
        ind_dir = np.argmax(full_like)
        cpmax_dir = CP.flatten()[ind_dir]
        cnmax_dir = CN.flatten()[ind_dir]
        fmax_dir = F.flatten()[ind_dir]

    if (mx > 1 and maj and PLOTLIKE):
        f, (ax1,ax2) = pl.subplots(2, figsize=(5, 9))
        pl.suptitle("N_grid = " + str(Ngrid))
        cf1 = ax1.contourf(np.log10(CP[:,:,0]), np.log10(CN[:,:,0]), full_like[:,:,0] - np.max(full_like),np.linspace(-1000,1,101))
        ax1.plot(np.log10(cpmax_maj_minus), np.log10(cnmax_maj_minus), 'gs')
        ax1.set_title("Negative")
        
        cf2 = ax2.contourf(np.log10(CP[:,:,-1]), np.log10(CN[:,:,-1]), full_like[:,:,-1] - np.max(full_like),np.linspace(-1000,1,101))
        ax2.plot(np.log10(cpmax_maj_plus), np.log10(cnmax_maj_plus), 'gs')
        ax2.set_title("Positive")
        pl.show()
    
    if (refine):
        #List of maximum likelihoods over each refinement step
        majlike_minus = np.zeros(Ndvals)
        majlike_plus = np.zeros(Ndvals)
        dirlike = np.zeros(Ndvals)
        
        #Refine for Majorana- and Dirac-like couplings
        #Based on current max-like values
        if (maj): #Majorana assumption
            cp1 = cpmax_maj_minus
            cn1 = cnmax_maj_minus
            f1 = -1.0 
            for i in range(Ndvals):
                (majlike_minus[i], cp1, cn1, f1)  = CalcLike_refine(mx, expts, Ngrid, cp1, cn1, f1, dvals[i], maj=True)
            
            cp1 = cpmax_maj_plus
            cn1 = cnmax_maj_plus
            f1 = 1.0
            for i in range(Ndvals):
                (majlike_plus[i], cp1, cn1, f1) = CalcLike_refine(mx, expts, Ngrid, cp1, cn1, f1, dvals[i], maj=True)
            
            maxlike = np.maximum(np.max(majlike_minus), np.max(majlike_plus))
        else: #Dirac assumption
            cp1 = cpmax_dir
            cn1 = cnmax_dir
            f1 = fmax_dir
            for i in range(Ndvals):
                (dirlike[i], cp1, cn1, f1) = CalcLike_refine(mx, expts, Ngrid, cp1, cn1, f1, dvals[i], maj=False)  
                
            maxlike = np.max(dirlike)
        
    else:
        maxlike = np.max(full_like)
    
    return maxlike
        
#-----------------
def CalcLike_refine(mx, expts, Ngrid, cp0, cn0, f0, delta, maj):
    
    #Range of values in f to grid scan
    deltaf = delta/4.0
    
    #Make sure we don't stray too far...
    if ((np.log10(cp0)-delta) < -11):
        cp0 = 1e-11*(10**delta)
    if ((np.log10(cn0)-delta) < -11):
        cn0 = 1e-11*(10**delta)
        
    #New grid based on max-like values
    cp_list = np.logspace(np.log10(cp0)-delta, np.log10(cp0)+delta, Ngrid)
    cn_list = np.logspace(np.log10(cn0)-delta, np.log10(cn0)+delta, Ngrid)

    N_expt = len(expts)
    if (maj):
        #Just sample case of f = +- 1
        f_list = np.asarray([f0])
        Nfvals = 1
    else:
        #Sample near max-like value 
        #(but not outside f = [-1, 1])
        fmin = f0-deltaf
        fmin = np.clip(fmin, -1.0, 1.0)
        
        fmax = f0+deltaf
        fmax = np.clip(fmax, -1.0, 1.0)
        
        f_list = np.linspace(fmin,fmax, Ngrid)
        Nfvals = Ngrid

    #Same procedure as in CalcLike_grid
    (CP, CN, F) = np.meshgrid(cp_list, cn_list, f_list, indexing='ij')
    
    full_like = 0.0
    for expt in expts:
        No = len(expt.events)
        like = 0.0
        A = np.zeros((expt.N_iso, Ngrid, Ngrid, Nfvals))
    
        for i in range(expt.N_iso):
            A[i, :, :, :] = 2.0*((CP*expt.N_p[i] + CN*expt.N_n[i])**2\
                     + 2*CP*CN*(F-1.0)*expt.N_n[i]*expt.N_p[i])

        A = np.clip(A, 1e-50, 1e50)

        if (expt.N_iso == 1):
            like = -A[0,:,:,:]*expt.Ne_list
            like += expt.eventlike + No*np.log(A[0,:,:,:])
        else:
            Ne_tot = np.einsum('ijkl,i', A, expt.Ne_list)
            like = -Ne_tot + No*np.log(Ne_tot)
            Nek = np.einsum('ijkl,i->ijkl', A, expt.Ne_list)/Ne_tot
            like += expt.eventlike_interp.ev(Nek[0,:,:,:],Nek[1, :,:,:])

        full_like += like

    #Get best fit for Majorana-like case
    ind_maj_minus = np.argmax(full_like[:,:,0].flatten())
    cpmax_maj_minus = CP[:,:,0].flatten()[ind_maj_minus]
    cnmax_maj_minus = CN[:,:,0].flatten()[ind_maj_minus]
    
    ind_maj_plus = np.argmax(full_like[:,:,-1].flatten())
    cpmax_maj_plus = CP[:,:,-1].flatten()[ind_maj_plus]
    cnmax_maj_plus = CN[:,:,-1].flatten()[ind_maj_plus]
   
    #Get best fit for Dirac-like case
    ind_dir = np.argmax(full_like)
    cpmax_dir = CP.flatten()[ind_dir]
    cnmax_dir = CN.flatten()[ind_dir]
    fmax_dir = F.flatten()[ind_dir]

    if (mx > 1 and maj and f0 < 0 and PLOTLIKE):
        f, (ax1,ax2) = pl.subplots(2, figsize=(5, 9))
        pl.suptitle("Refined")
        cf1 = ax1.contourf(np.log10(CP[:,:,0]), np.log10(CN[:,:,0]), full_like[:,:,0] - np.max(full_like),np.linspace(-1000,1,101))
        ax1.plot(np.log10(cpmax_maj_minus), np.log10(cnmax_maj_minus), 'gs')
        ax1.set_title("Negative")
        
        ax1.set_xlim(np.log10(cp0)-delta,np.log10(cp0)+delta)
        ax1.set_ylim(np.log10(cn0)-delta,np.log10(cn0)+delta)

        pl.show()
        
    # Return the maximum likelihood at this refinement step,
    # as well as best fit point
    return np.max(full_like), cpmax_dir, cnmax_dir, fmax_dir

#Calculate the significance given L0 and L1
def CalcSignificance(L0, L1):
    q = 2.0*(L1 - L0)
    #q is half-chi-squared distributed (because |f|=1 is on the boundary)
    pval = 0.5*(1-chi2.cdf(q,1))
    #Calculate a one-sided Gaussian significance
    sig = norm.ppf(1-pval)
    return sig
