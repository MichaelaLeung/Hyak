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
    cirrus = clouds/_cirrus_hitran2012_5000_20000cm_toa.rad
    cirrus_wl = cirrus.lam
    cirrus_flux = cirrus.pflux/cirrus_sflux
    strato = clouds/_strato_hitran2012_5000_20000cm_toa.rad
    strato_wl = strato.lam
    strato_flux = strato.pflux/strato_sflux
    avg_wl = (cirrus_wl[:len(strato_wl)] + strato_wl)/2
    avg_flux = (cirrus_flux[:len(strato_flux)] + strato_flux)/2
    fig, ax = plt.subplots(figsize = (30, 10))
    ax.plot(avg_wl, avg_flux)
    fig.savefig("avg_clougs.png", bbox_inches = 'tight')
    
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
        








