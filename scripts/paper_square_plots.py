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

def plotting(lamin, lamax, atmos, title):
    matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})
    matplotlib.rcParams['font.size'] = 25.0
    matplotlib.rc('text', usetex=False)
    plt.switch_backend('agg')
    fig_name = int(100*(float(lamin) + float(lamax))/2)
    wl, flux = cloud_weight(lamin, lamax, 0.01)
    if atmos == 0: # zero = ocean loss
        wl4, flux4 = ocean_loss(lamin, lamax)
        fig, ax = plt.subplots(figsize = (10,10))
        ax.plot(wl, flux, label = "1 bar Earth-Like")
        ax.plot(wl4, flux4, label = "10 bar Ocean Loss")
        ax.set_title(title)
        ax.set_ylabel("Reflectance")
        ax.set_xlabel("Wavelength ($\mu$m)")
        if lamin == 0.61:
            ax.legend()
        fig.savefig("/gscratch/vsm/mwjl/projects/high_res/plots/" + str(fig_name) +  "new_CIA_clouds.png", bbox_inches = "tight")
    else:
	wl_ocean, flux_ocean = cloud_weight_highw(lamin, lamax, 0.01)
        fig, ax = plt.subplots(figsize = (10,10))
        ax.plot(wl, avg_flux, label = "1 bar Earth-Like")
        ax.plot(wl_ocean, avg_flux2, label = "10 bar Ocean Outgassing")
        ax.set_title(title)
        ax.set_ylabel("Reflectance")
        ax.set_xlabel("Wavelength ($\mu$m)")
        if lamin == 0.61:
            ax.legend()
        fig.savefig("/gscratch/vsm/mwjl/projects/high_res/plots/" + str(fig_name) +  "new_CIA_ocean_clouds.png", bbox_inches = "t$


if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="nor_cld",
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
        plotting(0.67,0.71,0,"Oxygen B band (0.69) Ocean Loss")
        plotting(0.74,0.78,0,"Oxygen A band (0.76) Ocean Loss")
        plotting(1.25,1.29,0,"1.27 Ocean Loss")
        plotting(0.61,0.645,1,"Gamma band (0.63) Ocean Outgassing")
        plotting(0.67,0.71,1,"Oxygen B band (0.69) Ocean Outgassing")
        plotting(0.74,0.78,1,"Oxygen A band (0.76) Ocean Outgassing")
        plotting(1.25,1.29,1,"1.27 Ocean Outgassing")
    else:
	plotting(0.61,0.645,1,"Gamma band (0.63) Ocean Outgassing")





