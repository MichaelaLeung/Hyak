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
    sim.open_outputs()
    wl = data.lam
    flux = data.pflux
    sflux = data.sflux

    adj_flux = flux/sflux * ((sim.smartin.radius / sim.smartin.r_AU) **2 )

    fig, ax = plt.subplots(figsize = (30, 10))
    ax.plot(wl, adj_flux)
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$ m)")
    ax.set_title(label)
    fig.savefig(str(atmos) + ".png", bbox_inches = 'tight')

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







