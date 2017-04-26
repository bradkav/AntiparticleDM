from DMUtils import *
from Experiment import *
import sys
import matplotlib.pyplot as pl
from scipy.stats import chi2, norm
import CalcParamPoint as CPP

ensemble = sys.argv[1]
m0 = float(sys.argv[2])
index = int(sys.argv[3])


if (len(sys.argv) > 4):
    output_folder = sys.argv[4]
else:
    output_folder = "results/"

target_sigma = 1e-46
if (m0 > 90):
    target_sigma = 1e-45

#Recalculate to give the correct numbers!
l = CPP.GeneratePoint_ind(index)
sig0 = (1.973e-14*1.973e-14)*4.0*(reduced_m(1.0, m0))**2.0/np.pi
sig = sig0*0.5*((l[0]*54.0 + l[1]*77.0)**2.0 + (l[2]*54.0 + l[3]*77.0)**2.0)
sig_p = sig/(54+77)**2
l *= np.sqrt(target_sigma/sig_p)


print l
#Should be 312.544 events!

#Sampling paramaters
loglsq_min = -11
loglsq_max = -7

Nvals = 50

#----Functions----

def CalcLike_refine(mx, Nvals, cp0, cn0, f0, maj=False):
    cp_list = np.logspace(np.log10(cp0)-0.5, np.log10(cp0)+0.5, Nvals)
    if (maj):
        cn_list = np.asarray([1e-15,])
        Nn = 1
    else:
        cn_list = np.logspace(np.log10(cn0)-0.5, np.log10(cn0)+0.5, Nvals)
        Nn = Nvals
        
    fmin = f0-0.25
    if (fmin < -1.0):
        fmin = -1.0
        fmax = -0.5
    else:
        fmax = f0+0.25
        if (fmax > 1.0):
            fmax = 1.0
            fmin = 0.5

    f_list = np.linspace(fmin,fmax, Nvals)
    #f_list = np.linspace(-1.0, 1.0, Nvals)
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


print " Loading experiments for ensemble", ensemble, "..."

if (ensemble == "A"):
    exptlist = ["Xenon2", "Argon", "Silicon"]
elif (ensemble == "B"):
    exptlist = ["Xenon2", "Argon", "Germanium"]
elif (ensemble == "C"):
    exptlist = ["Xenon2", "Argon", "CaWO4"]
elif (ensemble == "D"):
    exptlist = ["Xenon2", "Argon", "Germanium_half","CaWO4_half"]


#exptlist = ["Xenon2", "Argon", "Silicon"]
N_expt = len(exptlist)
expts = [ Experiment(exptlist[i] + ".txt") for i in range(N_expt)]

#print " Generating events..."


print " Calculating likelihoods..."
Nmvals = 20
mlist = np.logspace(np.log10(25), np.log10(100), Nmvals)
#print mlist
likelist_maj = np.zeros((Nmvals))
likelist_dir = np.zeros((Nmvals))

#Nlist = np.array((4, 25, 50, 100))

#100 is the correct number, but 50 is quick!

Nsamps = 50

sigvals = np.zeros(Nsamps)

for k in range(Nsamps):
    for expt in expts:
        expt.GenerateEvents(m0, l)
    for i, mi in enumerate(mlist):
        #print i+1, mi
        #for j in range(len(Nlist)):
        likelist_maj[i],likelist_dir[i] = CalcLike_grid(mi, 100, refine=True)
    L0 = np.nanmax(likelist_maj)
    L1 = np.nanmax(likelist_dir)
    deltaL = L1 - L0
    pval = 1-chi2.cdf(2*deltaL,1)
    #We are using a 2-sided convention for the significance, Z  
    sig = norm.ppf(1-pval/2.0)
    sigvals[k] = sig
    print " Sample", str(k+1), " - Discrimination significance (N_grid = 64 - with refinement):", sig, "sigma"

np.savetxt(output_folder + "Results_p" + str(index)+".txt",sigvals, \
    header="Ensemble "+ ensemble + ", m = " + str(m0) + ", lambda = "+str(l))


print sigvals
print np.median(sigvals)
#pl.figure()
#pl.hist(sigvals, np.linspace(0,6,12))
#pl.show()

