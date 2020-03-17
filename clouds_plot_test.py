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

from rv_plots import ocean_outgassing

def plotting(lamin, lamax, res):
    wl, flux = cloud_weight(lamin, lamax, res)
    wl_oo, flux_oo = cloud_weight_highw(lamin, lamax, res)
    wl_oo_clouds, flux_oo_clouds = ocean_outgassing(lamin, lamax, res)
    fig, ax = plt.subplots(figsize=(10,10))
    ax.plot(wl,flux)
    ax.plot(wl_oo, flux_oo)
    ax.plot(wl_oo_clouds, flux_oo_clouds_)
    ax.legend()

if __name__ == '__main__':

    import platform
    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="nor_plt",
                               subname="submit.csh",
                               workdir = "",
                               nodes = 1,
                               mem = "500G",
                               walltime = "108:00:00",
                               ntasks = 28,
                               account = "vsm",
                               submit = True,
                               rm_after_submit = True)
    elif platform.node().startswith("n"):
        # On a mox compute node: ready to run
        plotting(0.76, 0.78. 0.01)
    else:
        plotting(0.76, 0.78. 0.01)

