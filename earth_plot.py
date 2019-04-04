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

def earth(lamin, lamax):
    res = 1/(10*lamin)

    
    sim = smart.interface.Smart(tag = 'earth')
    infile = "earth_avg.pt"
    sim.smartin.alb_file = "composite1_txt.txt"
    
    sim.set_run_in_place() 
    sim.set_executables_automatically()

    sim.smartin.sza = 57
    sim.load_atmosphere_from_pt(infile, addn2 = False)


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

def earth_noO4(lamin, lamax):
    res = 1/(10*lamin)

    
    sim = smart.interface.Smart(tag = 'earth_noO4')
    infile = "earth_avg.pt"
    sim.smartin.alb_file = "composite1_txt.txt"
    
    sim.set_run_in_place() 
    sim.set_executables_automatically()

    sim.smartin.sza = 57
    sim.load_atmosphere_from_pt(infile, addn2 = False)


    sim.smartin.FWHM = res
    sim.smartin.sample_res = res

    sim.smartin.minwn = 1e4/lamax
    sim.smartin.maxwn = 1e4/lamin 

    sim.lblin.minwn = 1e4/lamax
    sim.lblin.maxwn = 1e4/lamin

    o2 = sim.atmosphere.gases[2]
    o2.cia_file = None


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


def plotting(lamin, lamax):
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

    wl, flux = earth(lamin, lamax)
    wl2, flux2 = earth_noO4(lamin, lamax)
    fig, ax = plt.subplots(figsize = (10,10))
    ax.plot(wl, flux, label = "Earth")
    ax.plot(wl2, flux2, label = "Earth, no O2O2")
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$ m)")
    ax.legend()
    ax.set_title("Earth O2O2 test")
    fig.savefig("earth_o2o2test.png", bbox_inches = 'tight')

if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="earth",
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
        plotting(0.5,0.65)
    else:
        # Presumably, on a regular computer: ready to run
        plotting(0.55, 0.59)







