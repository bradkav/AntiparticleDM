import sys
import CalcParamPoint as CPP
from CalcLikelihood import *
from WIMpy.Experiment import Experiment
from WIMpy.MyParser import ReadParam
from scipy.interpolate import interp1d
from WIMpy.DMUtils import *
from scipy.integrate import quad

print " "
print " **********************************"
print " *      CalcDiscrimination.py     *"
print " **********************************"
print " Code for calculating significance of discriminating between Dirac and Majorana DM..."
print " Version 1.0.2 - BJK 06/07/2017"
print " "


def CalcDiscrim(ensemble, m0, f, r_np, outfile=None, exposure=-1):
    """
    CalcDiscrim(ensemble, m0, f, r_np, outfile, exposure=-1)
    
    Arguments:
        ensemble    - experimental ensemble to use: "A", "B", "C" or "D"
        m0          - input DM mass in GeV
        f           - interference term (Eq. 6). Note f = [-1, 1]
        r_np        - ratio of neutron to proton couplings lambda_n/lambda_p (Eq. 4 and 5)
                      NB: r_np > 0.
        outfile     - output file for the list of significances
        exposure    - exposure (in kg years) for experiments other than Xe and Ar. 
                      Set to -1 to ignore and just use input from files in DDexpt folder.
    
    """
    
    save_to_file = False
    if (outfile):
        save_to_file = True

    #Number of mock data sets to generate
    Nsamps = int(ReadParam("params.txt", "N_samples"))
    if (Nsamps == -1):
        Nsamps = 100

    #Select the correct experimental ensemble
    print " Loading experiments for ensemble", ensemble, "..."
    if (ensemble == "A"):
        exptlist = ["Xenon", "Argon", "Silicon"]
    elif (ensemble == "A_full"):
        exptlist = ["Xenon_full", "Argon", "Silicon_full"]
    elif (ensemble == "B"):
        exptlist = ["Xenon", "Argon", "Germanium"]
    elif (ensemble == "C"):
        exptlist = ["Xenon", "Argon", "CaWO4"]
    elif (ensemble == "D"):
        exptlist = ["Xenon", "Argon", "Germanium_half","CaWO4_half"]
    elif (ensemble == "D_full"):
        exptlist = ["Xenon_full", "Argon", "Germanium_half", "CaWO4_half"]

    exptdir = "DDexpt/"

    #Initialise experiments
    N_expt = len(exptlist)
    expts = [ Experiment(exptdir + exptlist[i] + ".txt") for i in range(N_expt)]

    #Generate couplings from values of r_np and f
    l = CPP.GeneratePoint(r_np, f)

    #This is where we calculate the normalisation of the couplings
    #Aim for the same number of events in Xenon as:
    # sigma = 1e-46 and m = 50
    sig0 = expts[0].sig_eff(50, l)
    l *= np.sqrt(1e-46/sig0)
    Ntarget = expts[0].CalcNevents(50, l)
    #Rescale to give correct event numbers for mass m0
    l *= np.sqrt(Ntarget/expts[0].CalcNevents(m0, l))

    #l is in the format [lpD, lnD, lpDb, lnDb]
    print " "
    print " DM mass [GeV]:", m0
    print " lambda [GeV^-2]:", l
    print " f =", f
    print " c_n/c_p =", r_np
    print " "
    for i in range(N_expt):
        print " Ne(" + exptlist[i] + "): ", expts[i].CalcNevents(m0,l), "; sig_p =", expts[i].sig_eff(m0, l)
    print " "

    #Set the exposure the experiments other than Xe and Ar (if required)
    if (exposure > 0):                           
        #Correction for number of experiments
        corr = 1.0/(N_expt - 2.0)
        for i in range(2, N_expt):
            expts[i].exposure = corr*exposure*365*0.7

        print " exposure (for other expts) = "+str(expts[-1].exposure/365.0) + " kg yr"
        print " "

    #Pre-tabulate some of the likelihood stuff:
    for expt in expts:
        expt.TabulateAll(m0)

    #Check against current limits - DarkSide-50
    DSdata = np.loadtxt("DS50limits_1510.00702.txt")
    DSlim = interp1d(DSdata[:,0], DSdata[:,1])

    exptArgon = Experiment(exptdir + "Argon.txt")
    if (exptArgon.sig_eff(m0, l) > DSlim(m0)):
        print " Current limit exceeded (DS-50)..."
        if (save_to_file):
            np.savetxt(outfile, sigvals, header="Ensemble "+ ensemble + ", m = " + str(m0) + ", lambda = "+str(l) + "; CURRENT LIMITS EXCEEDED (DS-50)")
        sys.exit()

    #Check against current limits - SuperCDMS
    CDMSdata = np.loadtxt("SuperCDMSlimits_1504.05871.txt")
    CDMSlim = interp1d(CDMSdata[:,0], CDMSdata[:,1])
    exptGermanium = Experiment(exptdir + "Germanium.txt")
    if (exptGermanium.sig_eff(m0, l) > CDMSlim(m0)):
        print " Current limit exceeded (SuperCDMS)..."
        if (save_to_file):
            np.savetxt(outfile, sigvals, header="Ensemble "+ ensemble + ", m = " + str(m0) + ", lambda = "+str(l) + "; CURRENT LIMITS EXCEEDED (SuperCDMS)")
        sys.exit()

    #List of masses to calculate for
    Nmvals = 25
    delta_m = np.log10(m0-10.0) #Size of range of masses depends on underlying DM mass
    m_min = m0/delta_m
    m_max = m0*delta_m
    if (m_min < 20.0):
        m_min = 20.0

    if (m0 >= 500):
        m_min = m0/10.0
        m_max = m0*10.0

    mlist = np.logspace(np.log10(m_min), np.log10(m_max), Nmvals)

    likelist_maj = np.zeros((Nmvals))
    likelist_dir = np.zeros((Nmvals))


    print (" Generating " + str(Nsamps) + " samples...")

    #Discrimination significance
    sigvals = np.zeros(Nsamps)


    #Generate Nsamps mock data sets
    for k in range(Nsamps):
        #Generate sample of events
        for expt in expts:
            expt.GenerateEvents(m0, l)
        #Calculate likelihood on a grid for range of mass values
        for i, mi in enumerate(mlist):
            #Majorana hypothesis
            likelist_maj[i] = CalcLike_grid(mi, expts, 200, maj=True, refine=True)
            #Dirac hypothesis
            likelist_dir[i] = CalcLike_grid(mi, expts, 50, maj=False, refine=True)
        sig = CalcSignificance(np.nanmax(likelist_maj), np.nanmax(likelist_dir))
        sigvals[k] = sig
        print "   Sample", str(k+1), " - Discrimination significance:", sig, "sigma"

    print " Median significance:", np.median(sigvals)

    if (exposure > 0):
        expstr = "exposure (of 3rd experiment) = " + str(exposure) + "kg years (before efficiency cut)"
    else:
        expstr = " "

    #Header information for the file
    hdr = '\n'.join(["Ensemble "+ ensemble + ", m = " + str(m0) + ", c_n/c_p = " + str(r_np) + ", f = " + str(f),\
                    "lambda = "+str(l),\
                     expstr])
                    
    #Output to file
    if (save_to_file):
        np.savetxt(outfile ,sigvals, header=hdr)
