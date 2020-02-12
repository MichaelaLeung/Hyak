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

def long():
    lamin = 0.4
    lamax = 2.0
    
    fig, ax = plt.subplots(figsize = (30,10))

    w_temp = smart.readsmart.Rad("highw_5000_25000cm_toa.rad")
    ocean_wl = w_temp.lam
    ocean_flux = w_temp.pflux
    ocean_sflux = w_temp.sflux
    ocean_flux = ocean_flux / ocean_sflux

    w_cir_temp = smart.readsmart.Rad("highw_cirrus_5000_25000cm_toa.rad")
    ocean_wl2 = w_cir_temp.lam
    ocean_flux2 = w_cir_temp.pflux
    ocean_sflux2 = w_cor_temp.sflux
    ocean_flux2 = ocean_flux2 / ocean_sflux2

    w_str_temp = smart.readsmart.Rad("highw_strato_5000_25000cm_toa.rad")
    ocean_wl3 = w_str_temp.lam
    ocean_flux3 = w_str_temp.pflux
    ocean_sflux3 = w_str_temp.sflux
    ocean_flux3 = ocean_flux3 / ocean_sflux3
    
    m, m_clouds = smart.utils.get_common_masks(ocean_wl, ocean_wl2)
    avg_flux = (0.5*ocean_flux[m_clouds]+0.25*ocean_flux2[m_clouds]+0.25*ocean_flux3[m_clouds])
    ax.plot(wl, avg_flux)
    ax.set_title("Simulated 10 bar oxygen ocean planet orbiting Proxima Centauri")
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$m)")


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
        long()
    else:
        long()
