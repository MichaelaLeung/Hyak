import numpy as np
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from astropy.io import fits
import smart
import sys, os
import datetime
matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['legend.loc'] = 'lower center'
import random
import math
import imageio
import platform 
import scipy
from scipy import interpolate
print('import finished') 
import scipy.integrate as integrate

from clouds_2020 import cloud_weight_highw
from clouds_2020 import cloud_weight

from rv_plots import ocean_loss
from rv_plots import run_earth

def plotting(atmos):
    lamin = 0.5
    lamax = 2.0
    res = 0.01
    
    import platform
    matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})
    matplotlib.rcParams['font.size'] = 25.0
    matplotlib.rc('text', usetex=False)
    plt.switch_backend('agg')
        
    fig, ax = plt.subplots(figsize = (30, 10))
    if atmos == 'prox': 
        wl, flux = cloud_weight(lamin, lamax, res)
        label = "Simulated Earth-like planet orbiting Proxima Centauri"
    elif atmos == 'highd':
        wl, flux = ocean_loss(lamin, lamax,res)
        label = "Simulated clear sky post ocean-loss planet orbiting Proxima Centauri"
    elif atmos == 'highw':
        wl, flux = cloud_weight_highw(lamin, lamax, res)
        label = "Simulated cloudy ocean outgassing planet orbiting Proxima Centauri"

    ax.plot(wl[:len(flux)], flux[:len(wl)])
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$m)")
    ax2 = ax.twinx()
    ax2.set_ylabel("Planet-to-star contrast ratio")        
    ax.set_title(label)
    ax.set_xlim(0.5,2)
    ax.axvspan(0.61, 0.65, alpha=0.5, color='0.85')
    ax.axvspan(0.67, 0.71, alpha=0.5, color='0.85')
    ax.axvspan(0.74, 0.78, alpha=0.5, color='0.85')
    ax.axvspan(1.25, 1.29, alpha=0.5, color='0.85')
    fig.savefig("/gscratch/vsm/mwjl/projects/high_res/plots/" + str(atmos) + "_newCIA.png", bbox_inches = 'tight')

if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="co_test",
                               subname="submit.csh",
                               workdir = "",
                               nodes = 1,
                               mem = "500G",
                               walltime = "24:00:00",
                               ntasks = 28,
                               account = "vsm",
                               submit = True,
                               rm_after_submit = True)
    elif platform.node().startswith("n"):
        # On a mox compute node: ready to run
     #   plotting("prox")
#        plotting("highd") 
#        plotting("highw")
#    else:
        # Presumably, on a regular computer: ready to run
        plotting("prox")      









