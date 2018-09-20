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

def smart_basic(res, lamin, lamax):
    res = 1/(10*lamin)
    sim = smart.interface.Smart(tag = "prox")
    sim2 = smart.interface.Smart(tag = "earth")
    infile = "10bar_O2_dry.pt_filtered.pt"
    earth_infile = "earth_avg.pt"
    
    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, "oxygen_outplot")

    try:
        os.mkdir(place)
    except OSError:
        pass

    sim.set_run_in_place(place)
        
    sim.set_executables_automatically()
    sim.load_atmosphere_from_pt(infile, addn2 = False)
    sim.set_planet_proxima_b()
    sim.set_star_proxima()
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

    flux = flux/sflux



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


    fig, ax = plt.subplots(figsize = (10,10))
    ax.plot(wl, flux)
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength (microns")
    fig_name = str(lamin) + "to" + str(lamax)
    fig.savefig("higho2" + fig_name +  ".png", bbox_inches = "tight")    
 

if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="manyplt",
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
        smart_basic(0.01, 0.75, 0.77)
        smart_basic(0.01, 1.26, 1.28)

    else:
        # Presumably, on a regular computer: ready to run
        smart_basic(0.01, 0.75, 0.77)
        smart_basic(0.01, 1.26, 1.28)
        smart_basic(0.01, 0.61, 0.65)
        smart_basic(0.01, 0.66, 0.70)








