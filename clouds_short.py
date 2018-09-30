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

def plotting():
    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, "clouds")

    try:
        os.mkdir(place)
    except OSError:
        pass

    sim = smart.interface.Smart(tag = "earth")
    sim.smartin.alb_file = "composite1_txt.txt"
    infile = "earth_avg.pt"

    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, "clouds")

    try:
        os.mkdir(place)
    except OSError:
        pass

        
    sim.set_run_in_place(place) 
    sim.set_executables_automatically()
    


    sim.load_atmosphere_from_pt(infile, addn2 = False)
    sim.set_planet_proxima_b()
        

    import platform
    
    if platform.system() == 'Darwin':
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
    cirrus = smart.readsmart.Rad("clouds/_cirrus_hitran2012_5000_20000cm_toa.rad")
    cirrus_wl = cirrus.lam
    cirrus_flux = cirrus.pflux
    cirrus_sflux = cirrus.sflux
    cirrus_flux = cirrus_flux/cirrus_sflux * ((sim.smartin.radius / sim.smartin.r_AU) **2 )
    strato = smart.readsmart.Rad("clouds/_strato_hitran2012_5000_20000cm_toa.rad")
    strato_wl = strato.lam
    strato_flux = strato.pflux
    strato_sflux = strato.sflux
    strato_flux = strato_flux/strato_sflux * ((sim.smartin.radius / sim.smartin.r_AU) **2 )
    avg_wl = (cirrus_wl[:len(strato_wl)] + strato_wl)/2
    avg_flux = (cirrus_flux[:len(strato_flux)] + strato_flux)/2
    fig, ax = plt.subplots(figsize = (30, 10))
    ax.plot(avg_wl, avg_flux)
    from matplotlib import rcParams
    matplotlib.rcParams.update({'font.size': 30})
    matplotlib.rc('xtick', labelsize=25) 
    matplotlib.rc('ytick', labelsize=25)
    matplotlib.rc('label', labelsize=25) 
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$ m)")
    ax.set_title("Earth")
    fig.savefig("avg_clouds.png", bbox_inches = 'tight')
    
if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="clouds",
                               subname="submit.csh",
                               workdir = "",
                               nodes = 1,
                               mem = "500G",
                               walltime = "5:00:00",
                               ntasks = 28,
                               account = "vsm",
                               submit = True,
                               rm_after_submit = True)
    elif platform.node().startswith("n"):
        # On a mox compute node: ready to run
        plotting()
    else:
        # Presumably, on a regular computer: ready to run
        cirrus_wl, cirrus_flux = clouds(10, 1, 0)
        strato_wl, strato_flux = clouds(10, 0, 1)
        avg_wl = (cirrus_wl + strato_wl)/2
        avg_flux = (cirrus_flux + strato_flux)/2
        fig, ax = plt.subplots(figsize = (30, 10))
        ax.plot(avg_wl, avg_flux)
        fig.savefig("avg_clougs_low.png", bbox_inches = 'tight')
        








