from DMUtils import *
from Experiment import *
import sys
#import matplotlib.pyplot as pl
from scipy.stats import chi2, norm


#Sampling paramaters
logc_min = -12
logc_max = -7

delta = 0.5


#----Functions----
        
def CalcLike_grid(mx, expts, Ngrid = 100, refine=False):
    
    N_expt = len(expts)

    #Tabulate rates for the specified mass
    for i in range(N_expt):
        expts[i].TabulateAll(mx)


    #Initial grid of parameters
    cp_list = np.logspace(logc_min, logc_max, Ngrid)
    cn_list = np.logspace(logc_min, logc_max, Ngrid)
    f_list = np.linspace(-1.0,1.0, Ngrid)

    (CP, CN, F) = np.meshgrid(cp_list, cn_list, f_list, indexing='ij')

    
    full_like = 0.0
    for expt in expts:
        No = len(expt.events)
        like = 0.0
        A = np.zeros((expt.N_iso, Ngrid, Ngrid, Ngrid))
    
        for i in range(expt.N_iso):
            A[i, :, :, :] = 2.0*((CP*expt.N_p[i] + CN*expt.N_n[i])**2\
                     + 2.0*CP*CN*(F-1)*expt.N_p[i]*expt.N_n[i])

        A = np.clip(A, 1e-50, 1e50)

        if (expt.N_iso == 1):
            like = -A[0,:,:,:]*expt.Ne_list
            like += expt.eventlike + No*np.log(A[0,:,:,:])
        else:
            like = -np.dot(A.T,expt.Ne_list).T
            like += np.sum(np.log(np.dot(A.T,expt.R_list).T), axis=0)

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

    if (mx > 1e30):
        f, (ax1,ax2) = pl.subplots(2, figsize=(5, 9))
        cf1 = ax1.contourf(np.log10(CP[:,:,0]), np.log10(CN[:,:,0]), full_like[:,:,0] - np.max(full_like[:,:,0]),np.linspace(-50,1,101))
        ax1.plot(np.log10(cpmax_maj_minus), np.log10(cnmax_maj_minus), 'gs')
        ax1.set_title("Negative")
        yvals = [np.log10(cnmax_maj_minus)-delta,np.log10(cnmax_maj_minus)+delta, np.log10(cnmax_maj_minus)+delta, np.log10(cnmax_maj_minus)-delta,np.log10(cnmax_maj_minus)-delta]
        xvals = [np.log10(cpmax_maj_minus)-delta,np.log10(cpmax_maj_minus)-delta, np.log10(cpmax_maj_minus)+delta, np.log10(cpmax_maj_minus)+delta,np.log10(cpmax_maj_minus)-delta]
        ax1.plot(xvals,yvals,'k-')
        #pl.colorbar(cf1)
        
        cf2 = ax2.contourf(np.log10(CP[:,:,-1]), np.log10(CN[:,:,-1]), full_like[:,:,-1] - np.max(full_like[:,:,-1]),np.linspace(-50,1,101))
        ax2.plot(np.log10(cpmax_maj_plus), np.log10(cnmax_maj_plus), 'gs')
        ax2.set_title("Positive")
        yvals = [np.log10(cnmax_maj_plus)-delta,np.log10(cnmax_maj_plus)+delta, np.log10(cnmax_maj_plus)+delta, np.log10(cnmax_maj_plus)-delta,np.log10(cnmax_maj_plus)-delta]
        xvals = [np.log10(cpmax_maj_plus)-delta,np.log10(cpmax_maj_plus)-delta, np.log10(cpmax_maj_plus)+delta, np.log10(cpmax_maj_plus)+delta,np.log10(cpmax_maj_plus)-delta]
        ax2.plot(xvals,yvals,'k-')
        #pl.colorbar(cf2)
        pl.show()
        
    if (refine):
        #Refine for Majorana- and Dirac-like couplings
        #Based on current max-like values
        reflike_maj_minus = CalcLike_refine(mx, expts, Ngrid+1, cpmax_maj_minus, cnmax_maj_minus, -1.0, maj=True)
        reflike_maj_plus = CalcLike_refine(mx, expts, Ngrid+1, cpmax_maj_plus, cnmax_maj_plus, 1.0, maj=True)
        reflike_dir = CalcLike_refine(mx, expts, Ngrid+1, cpmax_dir, cnmax_dir, fmax_dir, maj=False)
        #Need the +1 so that the original point is on the refined grid!
        
        reflike_maj = np.maximum(reflike_maj_minus, reflike_maj_plus)
        
        if (reflike_maj > reflike_dir):
            reflike_dir = reflike_maj
        return reflike_maj, reflike_dir
    else:
        L_maj = np.max(full_like[:,:,(0,-1)])
        L_dir = np.max(full_like)
    
        return L_maj, L_dir
        
#-----------------
def CalcLike_refine(mx, expts, Ngrid, cp0, cn0, f0, maj):
    
    #New grid based on max-like values
    cp_list = np.logspace(np.log10(cp0)-delta, np.log10(cp0)+delta, Ngrid)
    cn_list = np.logspace(np.log10(cn0)-delta, np.log10(cn0)+delta, Ngrid)

    N_expt = len(expts)

    if (maj):
        #Just sample case of f = +- 1
        f_list = np.asarray([f0])
        #f_list = np.asarray([-1.0, 1.0])
        Nfvals = 1
    else:
        #Sample near max-like value 
        #(but not outside f = [-1, 1])
        fmin = f0-0.25
        if (fmin < -1.0):
            fmin = -1.0
            fmax = -0.5
        else:
            fmax = f0+0.25
            if (fmax > 1.0):
                fmax = 1.0
                fmin = 0.5
        f_list = np.linspace(fmin,fmax, Ngrid)
        Nfvals = Ngrid

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
            like = -np.dot(A.T,expt.Ne_list).T
            #for j in range(No):
            like += np.sum(np.log(np.dot(A.T,expt.R_list).T), axis=0)
        full_like += like

    ind_plus = np.argmax(full_like[:,:,(-1)].flatten())
    ind_minus = np.argmax(full_like[:,:,(0)].flatten())

    ind = np.argmax(full_like)
    #print full_like.shape
    
    cpmax_maj_plus = CP[:,:,-1].flatten()[ind_plus]
    cnmax_maj_plus = CN[:,:,-1].flatten()[ind_plus]
    
    cpmax_maj_minus = CP[:,:,0].flatten()[ind_minus]
    cnmax_maj_minus = CN[:,:,0].flatten()[ind_minus]
    

    cpmax_dir = CP.flatten()[ind]
    cnmax_dir = CN.flatten()[ind]
    fmax_dir = F.flatten()[ind]
    #print cpmax_dir, cnmax_dir, fmax_dir
    #print "    Best-fit (Dir.):",fmax_dir

    #print "    delta-chi-sq:", -2*(np.max(full_like) - np.max(full_like[0,:,:]))
    #print "    Best-fit (Maj.):",cpmax_maj, cnmax_maj
    #print "    Best-fit (Dir.):",cpmax_dir, cnmax_dir
    #print " "
    if (mx > 1e30):
        f, (ax1,ax2) = pl.subplots(2, figsize=(5, 9))
        cf1 = ax1.contourf(np.log10(CP[:,:,0]), np.log10(CN[:,:,0]), full_like[:,:,0] - np.max(full_like[:,:,0]),np.linspace(-50,1,101))
        ax1.plot(np.log10(cpmax_maj_minus), np.log10(cnmax_maj_minus), 'gs')
        ax1.set_title("Negative")
        yvals = [np.log10(cnmax_maj_minus)-delta,np.log10(cnmax_maj_minus)+delta, np.log10(cnmax_maj_minus)+delta, np.log10(cnmax_maj_minus)-delta,np.log10(cnmax_maj_minus)-delta]
        xvals = [np.log10(cpmax_maj_minus)-delta,np.log10(cpmax_maj_minus)-delta, np.log10(cpmax_maj_minus)+delta, np.log10(cpmax_maj_minus)+delta,np.log10(cpmax_maj_minus)-delta]
        #ax1.plot(xvals,yvals,'k-')
        #pl.colorbar(cf1)
        
        cf2 = ax2.contourf(np.log10(CP[:,:,-1]), np.log10(CN[:,:,-1]), full_like[:,:,-1].T - np.max(full_like[:,:,-1]),np.linspace(-50,1,101))
        ax2.plot(np.log10(cnmax_maj_plus), np.log10(cpmax_maj_plus), 'gs')
        ax2.set_title("Positive")
        #pl.colorbar(cf2)
        pl.show()

    return np.max(full_like)
    
    #if (maj):
    #    return np.max(full_like)
    #else:
    #    return np.max(full_like)
    
def CalcSignificance(L0, L1):
    deltaL = L1 - L0
    pval = 1-chi2.cdf(2*deltaL,1)
    sig = norm.ppf(1-pval/2.0)
    return sig
