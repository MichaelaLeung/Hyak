import numpy as np
import smart
from matplotlib import pyplot as plt
import matplotlib as mpl
from matplotlib.collections import LineCollection
from astropy.io import fits 
import matplotlib
import sys, os
import datetime
matplotlib.rcParams['text.usetex'] = False
import random

def earth_like(lamin, lamax):
    res = 1/(10*lamin)

    sim = smart.interface.Smart(tag = "prox")
    infile = "profile_Earth_proxb_.pt_filtered"
    label = "Simulated Earth-like planet orbiting Proxima Centauri"
    sim.smartin.alb_file = "composite1_txt.txt"
    sim.set_planet_proxima_b()
    sim.load_atmosphere_from_pt(infile, addn2 = False)

    sim.lblin.par_file = '/gscratch/vsm/alinc/fixed_inputs/HITRAN2016.par'
    sim.lblin.hitran_tag = 'hitran2016'
    sim.lblin.fundamntl_file = '/gscratch/vsm/alinc/fixed_inputs/fundamntl2016.dat'
    sim.lblin.lblabc_exe = '/gscratch/vsm/alinc/exec/lblabc_2016'
    
    o2 = sim.atmosphere.gases[3]
    o2.cia_file = 'cia_adj_calc.cia'
    infile = "profile_Earth_proxb_.pt_filtered"
    label = "Earth-Like"
    sim.smartin.alb_file = "composite1_txt.txt"
    sim.set_planet_proxima_b()
    sim.set_star_proxima()

    sim.set_run_in_place() 
    sim.set_executables_automatically()

    sim.smartin.sza = 57

    sim.smartin.FWHM = res
    sim.smartin.sample_res = res

    sim.smartin.minwn = 1e4/lamax
    sim.smartin.maxwn = 1e4/lamin 

    sim.lblin.minwn = 1e4/lamax
    sim.lblin.maxwn = 1e4/lamin 


    sim.gen_lblscripts()
    sim.run_lblabc()
    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux
    adj_flux = flux/sflux
    return(wl, adj_flux)

def plotting(lamin, lamax, atmos, title):
    import platform
    if platform.system() == 'Jarvis':
        # On a Mac: usetex ok
        mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})
        mpl.rcParams['font.size'] = 25.0
        mpl.rc('text', usetex=True)
    elif platform.node().startswith("n"):
        # On hyak: usetex not ok, must change backend to 'agg'
        mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})
        mpl.rcParams['font.size'] = 25.0
        mpl.rc('text', usetex=False)
        plt.switch_backend('agg')
    fig_name = int(100*(float(lamin) + float(lamax))/2)
    if atmos == 0: # zero = ocean loss
        wl, flux = earth_like(lamin, lamax)
        wl2, flux2 = ocean_loss(lamin, lamax)
        fig, ax = plt.subplots(figsize = (10,10))
        ax.plot(wl, flux, label = "1 bar Earth-Like")
        ax.plot(wl2, flux2, label = "10 bar Ocean Loss")
        ax.set_title(title)
        ax.set_ylabel("Reflectance")
        ax.set_xlabel("Wavelength ($\mu$ m)")
        ax.legend()
        fig.savefig(str(fig_name) +  "2016_test.png", bbox_inches = "tight")
       
if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="smartplt",
                               subname="submit.csh",
                               workdir = "",
                               nodes = 1,
                               mem = "500G",
                               walltime = "10:00:00",
                               ntasks = 28,
                               account = "vsm",
                               submit = True,
                               rm_after_submit = True)
    elif platform.node().startswith("n"):
        # On a mox compute node: ready to run
        plotting(0.61,0.645,0,"Gamma band (0.63) Ocean Loss")
            else:
        plotting(0.61,0.645,1,"Gamma band (0.63) Ocean Outgassing")
