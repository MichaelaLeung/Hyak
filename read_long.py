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

def longplot(atmos):
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
        
    if atmos == "earth":
        infile = "earth_avg.pt"
    elif atmos == "prox":
        data = smart.readsmart.Rad("longplot/prox_hitran2012_5000_20000cm_toa.rad")
        label = "Self Consistent PCb (Earth like)"
    elif atmos == "highd":
        data = smart.readsmart.Rad("longplot/highd_hitran2012_5000_20000cm_toa.rad")
        label = "10 bar O2 PCb"
    elif atmos == "low":
        data = smart.readsmart.Rad("high_o2_noO4_hitran2012_11111_11627cm_toa.rad")
        label = "10 Bar O2 PCb"
    wl = data.lam
    flux = data.pflux
    sflux = data.sflux
    radius = 6850.0
    r_AU = 0.0485
    adj_flux = flux/sflux * ((radius / r_AU) **2 )

    refl = flux/sflux
    flux = adj_flux
    a = [1.25,1.25,1.27,1.27]
    b = [0, max(flux), max(flux), 0]
    c = [0.74,0.74,0.76,0.76]
    d = [0, max(flux), max(flux), 0]
    e = [0.61,0.61,0.65,0.65]
    f = [0, max(flux), max(flux), 0]
    g = [0.66,0.66,0.70,0.70]
    h = [0, max(flux), max(flux), 0]
    
    fig, ax = plt.subplots(figsize = (40, 10))
    ax.fill(a,b, '0.75')
    ax.fill(c,d, '0.75')
    ax.fill(e,f, '0.75')
    ax.fill(g,h, '0.75')
    ax.plot(wl, flux)
    ax2 = ax.twinx()
    ax2.plot(wl, refl, 'b')
    ax2.xaxis.set_visible(False)
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$ m)")
    ax.set_title(label)
    plt.xticks(np.arange(min(wl), max(wl)+.1, 0.1))
 #   ax2.xticks(np.arange(min(wl), max(wl)+.1, 0.1))
 #   ax3.xticks(np.arange(min(wl), max(wl)+.1, 0.1))
    ax.set_xlim(0.5,2)
    ax2.set_xlim(0.5,2)
 #   ax3.set_xlim(0.5,2)

    fig.savefig(str(atmos) + ".png", bbox_inches = 'tight')

if __name__ == '__main__':
    import platform
    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="longplot",
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
        longplot("prox")
        longplot("highd")
    else:
        # Presumably, on a regular computer: ready to run
        longplot("low")






