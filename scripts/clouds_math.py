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

high_res = 0.01
low_res = 1

def high_pass(wl, flux, lamin, lamax, wl_low, flux_low):
    long_flux = []
    temp_factor = math.ceil(len(wl) / len(wl_low))
    print(temp_factor)
    for i in flux_low:
        j = 0
        while j < temp_factor: 
            long_flux.append(i)
            j = j+1
    mixed = []
    i = 0
    while i < len(flux):
        temp = (flux[i] + long_flux[i]) / 2
        mixed.append(temp)
        i = i+1

    i = 0
    flattened = []
    while i < len(long_flux)- 25: 
        avg = np.mean(flux[i:i+25])
        j = 0
        while j < 25:
            flattened.append(avg)
            j = j+1
        i = i+25
    out = []
    i = 0
    while i < len(mixed[:-25]): 
        diff = abs(mixed[i] - flattened[i])
        out.append(diff)
        i = i+1
    fig, ax = plt.subplots(figsize = (10,10))
    ax.plt(wl, out)
    fig.save_fig('highpass'+str(lamin)+'.png')
    return(wl, out)

def flux_calc(lamin,lamax, type):
    print("flux calc", lamin)
    earth_wl, earth_flux = run_earth(lamin,lamax, high_res)
    earth_wl_low, earth_flux_low = run_earth(lamin, lamax,low_res)
    if type == 0:
         wl, flux = clouds_out(lamin, lamax,high_res)
         wl_low, flux_low = clouds_out(lamin, lamax, low_res)
    elif type == 1:
         wl, flux = ocean_loss(lamin,lamax,high_res)
         wl_low, flux_low = ocean_loss(lamin, lamax, low_res)
    elif type == 2:
         wl, flux = cloud_weight_highw(lamin,lamax, high_res)
         wl_low, flux_low = cloud_weight_highw(lamin, lamax,low_res)

    phases = np.linspace(0,2*np.pi,n_phase)
    inclination = np.pi/2
    phi_90 = np.pi/2
    sma = 7500000
    c = 299792.458

    fluxes = np.outer(earth_flux, np.ones(n_phase))
    temp = np.arccos(-np.sin(inclination)*np.cos(phases))
    phase_function = ((np.sin(temp)+(np.pi-temp)*(np.cos(temp)))/(np.pi))
    #season = 0.55 #0-2PI march max
    season = -0.55 #july min?
    rv_bary = (29.8 * np.sin((11.2/365.)*phases + season))
    #rv_orb = (np.sqrt((G*m_prox)/(sma_m)) * np.sin(inclination)/1000 *np.sin(phases))
    rv_orb = (sma*2*np.pi)/967680* np.sin(inclination) *np.sin(phases)
    rv_sys = -21.7 * np.ones_like(phases)
    #rv_bary = (29.8 * np.sin((11.2/365.)*phases))
    rv = rv_sys + rv_orb - rv_bary
    inclination = np.pi/2.9
    obs_wl = np.outer(wl,(1+rv/c))

    print(len(earth_wl))
    length = np.shape(earth_wl)
    length = length[0]
    out =np.empty((n_phase, length))

    i = 0
    while i < n_phase:
        high_pass_wl, high_pass_flux = high_pass(obs_wl[:len(flux),i],flux[:len(obs_wl)],lamin, lamax,wl_low, flux_low) 
        interp = scipy.interpolate.interp1d(high_pass_wl[:len(high_pass_flux)], high_pass_flux[:len(high_pass_wl)],fill_value = "extrapolate")
        temp = np.asarray(interp(earth_wl) * earth_flux)
        out[i,:] = (temp) 
        i = i+1
    fig, ax = plt.subplots(figsize = (10,10))
    ax.plt(wl, out)
    fig.save_fig('flux_calc'+str(lamin)+'.png')
    return(wl, out)

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
        print('job submitted') 
#        integ_calc(0.74, 0.78, 0)
#        read_integ()
        flux_calc(0.76, 0.78,0)
    else:
        flux_calc(0.60,0.70)

