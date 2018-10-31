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
    if atmos == "earth":
        infile = "earth_avg.pt"
    elif atmos == "prox":
        data = smart.readsmart.Rad("longplot/prox_hitran2012_5000_20000cm_toa.rad")
        label = "Self Consistent PCb"
    elif atmos == "highd":
        data = smart.readsmart.Rad("longplot/highd_hitran2012_5000_20000cm_toa.rad")
        label = "10 bar O2 PCb"
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
    
    fig, ax = plt.subplots(figsize = (30, 10))
    ax2 = ax.twiny()
    ax2.fill(a,b, '0.75')
    ax2.fill(c,d, '0.75')
    ax2.fill(e,f, '0.75')
    ax2.fill(g,h, '0.75')
    ax.plot(wl, flux)
    ax2.xaxis.set_visible(False)
    ax3 = ax.twinx()
    ax3.plot(wl, refl, 'r')
    ax2.xaxis.set_visible(False)
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$ m)")
    ax.set_title(label)
    ax.set_xlim(0.5,2)
    fig.savefig(str(atmos) + info + ".png", bbox_inches = 'tight')

if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="longplt",
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
        pass







