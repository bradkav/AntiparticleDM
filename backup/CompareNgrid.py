from DMUtils import *
from Experiment import *
import sys
import matplotlib.pyplot as pl
from scipy.stats import chi2, norm
import emcee
import CalcParamPoint as CPP

m0 = float(sys.argv[1])
index = int(sys.argv[2])

target_sigma = 1e-46
if (m0 > 90):
    target_sigma = 1e-45

#Recalculate to give the correct numbers!
l = CPP.GeneratePoint_ind(index)
sig0 = (1.973e-14*1.973e-14)*4.0*(reduced_m(1.0, m0))**2.0/np.pi
sig = sig0*0.5*((l[0]*54.0 + l[1]*77.0)**2.0 + (l[2]*54.0 + l[3]*77.0)**2.0)
sig_p = sig/(54+77)**2
l *= np.sqrt(target_sigma/sig_p)

l0 = l
#Should be 312.544 events!

print l0

#Sampling paramaters
loglsq_min = -11
loglsq_max = -7

Ngrid = 50


#----Functions----



print " Loading experiments..."
exptlist = ["Xenon2", "Argon", "Silicon"]
N_expt = len(exptlist)
expts = [ Experiment(exptlist[i] + ".txt") for i in range(N_expt)]

#expt1 = Experiment("Xenon.txt")
#print expt1.CalcNevents(m0,l0)
#expt1.GenerateEvents(m0, l0)
#expt1.PrintEvents()

print " Generating events..."
for i in range(N_expt):
    expts[i].GenerateEvents(m0, l0)
    #expts[i].TabulateRate(1e1, 1e3)


def lnprior(x):
    (lcp, lcn, f) = x
    if ((-11 < lcp < -5)and(-11 < lcn < -5)and(-1 < f < 1)):
        return 0.0
    return -np.inf

def CalcExptLike(expt, mx, l, T=1.0):

    #Poisson likelihood
    No = len(expt.events)
    Ne = expt.CalcNevents_tab(mx,l)
    #print Ne
    PL = -Ne + No*np.log(Ne) 

    #print expt, mx, l
    #print PL
    #Event-by-event likelihood
    for i in range(No):
        PL += np.log(expt.exposure*expt.dRdE(expt.events[i], mx, l)/Ne)

    if np.isnan(PL):
        return -1e30
    else:
        return PL*1.0/T

#Get a (sorted) chain with two columns - values and log-likelihoods
def ProcessChain(sampler):
    chain = np.column_stack((sampler.flatchain[:],sampler.flatlnprobability))
    chain = chain[chain[:,0].argsort()]
    return chain
        
"""
N1 = 50
N2 = 100
cp_list1 = np.logspace(loglsq_min, loglsq_max, N1)**2
cp_list2 = np.logspace(loglsq_min, loglsq_max, N2)**2



print len(range(1,101,1))
 


pl.figure()


pl.semilogy(np.arange(1,101, 2),cp_list1, '^')
pl.semilogy(np.arange(1, 101, 1),cp_list2, 's')

for i in range(50):
    pl.axhline(cp_list1[i], color='red')
for i in range(100):
    pl.axhline(cp_list2[i], color='blue')

pl.show()
#f_list = np.append(-1*np.logspace(0, -1, Nvals/2),np.logspace(-1, 0, Nvals/2))
"""
  
def CalcLike_refine(mx, Nvals, cp0, cn0, f0, maj=False):


    
    cp_list = np.logspace(np.log10(cp0)-1.0, np.log10(cp0)+1.0, Nvals)
    if (maj):
        cn_list = np.asarray([1e-15,])
        Nn = 1
    else:
        cn_list = np.logspace(np.log10(cn0)-1.0, np.log10(cn0)+1.0, Nvals)
        Nn = Nvals
    f_list = np.linspace(np.maximum(-1.0, f0-0.5),np.minimum(f0+0.5, 1.0), Nvals)
    #print f_list
    #print cp_list
    #print f_list
    #for i in range(N_expt):
    #    expts[i].TabulateAll(mx)

    (CPsq, CNsq, Fr) = np.meshgrid(cp_list**2, cn_list**2, f_list)
    
    Ident = np.ones((Nn, Nvals, Nvals))
    
    full_like = 0.0
    for expt in expts:
        No = len(expt.events)
        like = 0.0
        A = np.zeros((expt.N_iso, Nn, Nvals, Nvals))
    
        for i in range(expt.N_iso):
            A[i, :, :, :] = (CPsq*(Ident*expt.N_p[i] + Fr*expt.N_n[i])**2 + CNsq*expt.N_n[i]**2)
        
        #print A*expt.Ne_list
    
        #This definitely doesn't work!
        #like = -np.dot(A.T, expt.Ne_list)
        

        if (expt.N_iso == 1):
            like = -A[0,:,:,:]*expt.Ne_list
            like += expt.eventlike + No*np.log(A[0,:,:,:])
        else:
            print "No currently working for N_iso > 1"
            like = -np.sum(A, axis=0)*expt.Ne_list
            like += np.sum(np.log(np.dot(A.T,(expt.R_list).T)), axis=3)
        #print "4"
        full_like += like
    
    #print full_like.shape
    ind = np.argmax(full_like[0,:,:])
    ind2 = np.argmax(full_like)

    #print full_like.shape
    
    cpmax_maj = np.sqrt(CPsq.flatten()[ind])
    fmax_maj = Fr.flatten()[ind]

    cpmax_dir = np.sqrt(CPsq.flatten()[ind2])
    cnmax_dir = np.sqrt(CNsq.flatten()[ind2])
    fmax_dir = Fr.flatten()[ind2]
    
    #print "    delta-chi-sq:", -2*(np.max(full_like) - np.max(full_like[0,:,:]))
    #if (maj):
    #    print "    Best-fit (Maj.):", cpmax_maj, fmax_maj*cpmax_maj
    #print "    Best-fit (Dir.):",cpmax_dir, fmax_dir*cpmax_dir, cnmax_dir
    #print " "
    
    if (maj):
        return np.max(full_like[0,:,:])
    else:
        return np.max(full_like)

        
def CalcLike_grid(mx, Nvals = 100, refine=False):

    #cp_list = np.logspace(-9, -7, Nvals)**2
    cp_list = np.append(1e-15,np.logspace(loglsq_min, loglsq_max, Nvals))**2
    cn_list = np.append(1e-15,np.logspace(loglsq_min, loglsq_max, Nvals))**2
    #f_list = np.append(-1*np.logspace(0, -1, Nvals/2),np.logspace(-1, 0, Nvals/2))
    f_list = np.linspace(-1.0,1.0, Nvals)

    for i in range(N_expt):
        expts[i].TabulateAll(mx)

    (CPsq, CNsq, Fr) = np.meshgrid(cp_list, cn_list, f_list)
    
    Ident = np.ones((Nvals+1, Nvals+1, Nvals))
    
    full_like = 0.0
    for expt in expts:
        No = len(expt.events)
        like = 0.0
        A = np.zeros((expt.N_iso, Nvals+1, Nvals+1, Nvals))
    
        for i in range(expt.N_iso):
            A[i, :, :, :] = (CPsq*(Ident*expt.N_p[i] + Fr*expt.N_n[i])**2 + CNsq*expt.N_n[i]**2)
        
        #print A*expt.Ne_list
    
        #This definitely doesn't work!
        #like = -np.dot(A.T, expt.Ne_list)
        
        #print "   ",No, np.min(A*expt.Ne_list),np.max(A*expt.Ne_list)
        #print -like
        #print like
        #print "3"
        if (expt.N_iso == 1):
            like = -A[0,:,:,:]*expt.Ne_list
            like += expt.eventlike + No*np.log(A[0,:,:,:])
        else:
            print "No currently working for N_iso > 1"
            like = -np.sum(A, axis=0)*expt.Ne_list
            like += np.sum(np.log(np.dot(A.T,(expt.R_list).T)), axis=3)
        #print "4"
        full_like += like

    #print full_like.shape
    ind = np.argmax(full_like[0,:,:])
    ind2 = np.argmax(full_like)

    #print full_like.shape
    
    cpmax_maj = np.sqrt(CPsq.flatten()[ind])
    fmax_maj = Fr.flatten()[ind]

    cpmax_dir = np.sqrt(CPsq.flatten()[ind2])
    cnmax_dir = np.sqrt(CNsq.flatten()[ind2])
    fmax_dir = Fr.flatten()[ind2]
    
    #print "    delta-chi-sq:", -2*(np.max(full_like) - np.max(full_like[0,:,:]))
    #print "    Best-fit (Maj.):",cpmax_maj, fmax_maj*cpmax_maj
    #print "    Best-fit (Dir.):",cpmax_dir, fmax_dir*cpmax_dir, cnmax_dir
    #print " "

    if (refine):
        reflike_maj = CalcLike_refine(mx, 100, cpmax_maj, 1e-30, fmax_maj, maj=True)
        reflike_dir = CalcLike_refine(mx, 100, cpmax_dir, cnmax_dir, fmax_dir, maj=False)
        return reflike_maj, reflike_dir
    else:


        
        if (mx > 1e30):
            pl.figure()
            pl.contourf(np.log10(np.sqrt(cp_list)), f_list, full_like[0,0,:,:].T - np.max(full_like[0,0,:,:]),np.linspace(-50,1,101))
            pl.plot(np.log10(np.sqrt(CPsq.flatten()[ind])), Fr.flatten()[ind], 'gs')
            pl.title(r'$N_\mathrm{grid} = '+ str(Nvals)+'$')
            pl.colorbar()
            pl.show()
    

        #Use "MAP"?
    
        #Monte carlo each value of the mass...
    
        #print np.argmax(full_like[:,0,:,:])

        #print cn_list[0]
        #ind = np.argmax(full_like)
        #print (A.flatten())[ind]*expts[0].Ne_list
        #print ind
        #print CPsq.flatten()[ind], CNsq.flatten()[ind], Fr.flatten()[ind]
        L_maj = np.max(full_like[0,:,:])
        L_dir = np.max(full_like)
    
        return L_maj, L_dir


def lnprob(x):
    lp = lnprior(x)
    if not np.isfinite(lp):
        return -np.inf

    cp = 10**x[0]
    cn = 10**x[1]
    f = x[2]
    
    like_full = 0.0
    for expt in expts:
        like = 0.0
        A = np.zeros(expt.N_iso)
        for i in range(expt.N_iso):
            
            #A[i] = 2.0*(cp**2*expt.N_p[i]**2 + cn**2*expt.N_n[i]**2 + 2.0*cp*cn*f*expt.N_p[i]*expt.N_n[i])
            A[i] = 2.0*(cp**2*(expt.N_p[i] + f*expt.N_n[i])**2 + cn**2*expt.N_n[i]**2)
        like = -np.dot(A.T, expt.Ne_list)
        if (expt.N_iso == 1):
            like += expt.eventlike + np.log(A[0])*len(expt.events)
        else:
            like += np.sum(np.log(np.dot(A.T,(expt.R_list).T)))
        if (np.isnan(like)):
            return -np.inf
        like_full += like
    return lp + like_full/Tchain

def CalcExptLike_MC(mx):

    for i in range(N_expt):
        expts[i].TabulateAll(mx)

    
    pos = [(-6.1,-6.2, 0.1) + 0.25*np.random.randn(ndim) for i in range(nwalkers)]
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=())
    sampler.run_mcmc(pos, nsteps)
    #chain = ProcessChain(sampler)
    #print chain.shape
    #print np.max(sampler.flatlnprobability)
    return np.max(sampler.flatlnprobability)*Tchain
    

    
#E_list = np.linspace(expts[0].E_min, expts[0].E_max, 100)
#R_list = 0.0*E_list
#for i, Ei in enumerate(E_list):
#    R_list[i] = expts[0].dRdE_tab(Ei, m0, l0)

#pl.figure()
#pl.loglog(E_list, R_list)
#pl.ylim(1e-10,100)
#pl.show()

print " Calculating likelihoods..."
Nmvals = 20
mlist = np.logspace(np.log10(30), np.log10(80), Nmvals)
likelist = np.zeros((Nmvals, 5))

likelist_low = np.zeros((Nmvals))
likelist_high = np.zeros((Nmvals))
likelist_max = np.zeros((Nmvals))

likelist_maj = np.zeros((Nmvals))
likelist_dir = np.zeros((Nmvals))

likelist_maj2 = np.zeros((Nmvals))
likelist_dir2 = np.zeros((Nmvals))

likelist_maj3 = np.zeros((Nmvals))
likelist_dir3 = np.zeros((Nmvals))

likelist_maj4 = np.zeros((Nmvals))
likelist_dir4 = np.zeros((Nmvals))

#Nlist = np.array((4, 25, 50, 100))

#100 is the correct number, but 50 is quick!

N0 = 50
N1 = 100
N2 = 200
N3 = 101

for i, mi in enumerate(mlist):
    print i+1, mi
    #for j in range(len(Nlist)):
    likelist_maj[i],likelist_dir[i] = CalcLike_grid(mi, N0)
    print " N = 50          : ",likelist_maj[i],likelist_dir[i]
    
    likelist_maj2[i],likelist_dir2[i] = CalcLike_grid(mi, N1)
    print " N = 100         : ",likelist_maj2[i],likelist_dir2[i]
    
    likelist_maj3[i],likelist_dir3[i] = CalcLike_grid(mi, N2)
    print " N = 200         : ",likelist_maj3[i],likelist_dir3[i]
    
    likelist_maj4[i], likelist_dir4[i] = CalcLike_grid(mi, N0, refine=True)
    print " N = 50 (ref.)   : ",likelist_maj4[i],likelist_dir4[i]

    #likelist_maj4[i],likelist_dir4[i] = CalcLike_grid(mi, N3)
    #likelist_maj[i],likelist_dir[i] = (-1e30, -1e30)


    
L0 = np.nanmax(likelist_maj)
L1 = np.nanmax(likelist_dir)
deltaL = L1 - L0
pval = 1-chi2.cdf(2*deltaL,1)
sig = norm.ppf(1-pval/2.0)

L0_2 = np.nanmax(likelist_maj2)
L1_2 = np.nanmax(likelist_dir2)
deltaL = L1_2 - L0_2
pval = 1-chi2.cdf(2*deltaL,1)
sig2 = norm.ppf(1-pval/2.0)

L0_3 = np.nanmax(likelist_maj3)
L1_3 = np.nanmax(likelist_dir3)
deltaL = L1_3 - L0_3
pval = 1-chi2.cdf(2*deltaL,1)
sig3 = norm.ppf(1-pval/2.0)

L0_4 = np.nanmax(likelist_maj4)
L1_4 = np.nanmax(likelist_dir4)
deltaL = L1_4 - L0_4
pval = 1-chi2.cdf(2*deltaL,1)
sig4 = norm.ppf(1-pval/2.0)

#We are using a 2-sided convention for the significance, Z                         

print " Discrimination significance (N_grid = "+str(N0)+"):", sig, "sigma"
print " Discrimination significance (N_grid = "+str(N1)+"):", sig2, "sigma"
print " Discrimination significance (N_grid = "+str(N2)+"):", sig3, "sigma"
print " Discrimination significance (N_grid = "+str(N3)+"):", sig4, "sigma"


 
"""
pl.figure()
for j in range(len(Nlist)):
    pl.semilogx(mlist, -2*(likelist[:,j]-L0), label=r"$N_\mathrm{grid} = "+str(Nlist[j])+"$")

pl.legend(loc="best", frameon=False)
pl.ylim(-1, 100)
pl.axvline(50, linestyle='--', color='k')
pl.axhline(0, linestyle='--', color='k')
pl.show()
"""

pl.figure()
pl.semilogx(mlist, -2*(likelist_maj3-L1_3), 'b-',label=r"Majorana", linewidth=1.5)
pl.semilogx(mlist, -2*(likelist_dir3-L1_3), 'g-',label=r"Dirac", linewidth=1.5)
pl.semilogx(mlist, -2*(likelist_maj2-L1_3),  'b--',linewidth=1.5)
pl.semilogx(mlist, -2*(likelist_dir2-L1_3), 'g--',linewidth=1.5)
pl.semilogx(mlist, -2*(likelist_maj4-L1_4),  'b:',linewidth=1.5)
pl.semilogx(mlist, -2*(likelist_dir4-L1_4), 'g:',linewidth=1.5)

pl.semilogx(1e-30, 1e-30, 'k:',label=r"$N_\mathrm{grid}="+str(N0)+"$", linewidth=1.5)
pl.semilogx(1e-30, 1e-30, 'k--',label=r"$N_\mathrm{grid}="+str(N1)+"$", linewidth=1.5)
pl.semilogx(1e-30, 1e-30, 'k-',label=r"$N_\mathrm{grid}="+str(N2)+"$", linewidth=1.5)

pl.legend(loc="best", frameon=False)
pl.ylim(-1, 30)
pl.xlim(10, 1000)
pl.axvline(50, linestyle='--', color='k')
pl.axhline(0, linestyle='--', color='k')
pl.show()

#print likelist
#print likelist
#pl.figure()
#pl.semilogx(mlist, -2*(likelist - np.nanmax(likelist)))

#pl.axvline(50, linestyle='--', color='k')
#pl.show()

#pl.figure()
#pl.hist(expt1.events,100)
#pl.show()

