#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 18:36:23 2020

@author: mwl
"""

import smart 
import matplotlib.pyplot as plt 

def final_plot():
    loc_tag = '/gscratch/vsm/mwjl/projects/high_res/smart_output/'
    file1 = loc_tag + 'prox_cirrus0.01_5000_20000cm_toa.rad'
    file2 = loc_tag + 'prox_strato_5000_20000cm_toa.rad'
    file3 = loc_tag + 'prox0.01_5000_20000cm_toa.rad'
    data = smart.readsmart.read_rad(file1)
    wl = data.lam
    flux = data.pflux
    sflux = data.sflux
    
    data2 = smart.readsmart.read_rad(file2)
    wl2 = data2.lam
    flux2 = data2.pflux
    sflux2 = data2.sflux
    
    data3 = smart.readsmart.read_rad(file3)
    wl3 = data3.lam
    flux3 = data3.pflux
    sflux3 = data3.sflux
    
    m, m_clouds = smart.utils.get_common_masks(wl, wl2)
    avg_flux = (0.5*flux[m]+0.25*flux2[m]+0.25*flux3[m])
    avg_sflux = (0.5*sflux[m]+0.25*sflux2[m]+0.25*sflux3[m])
    
    lamin = 0.5
    lamax = 2.0
    res = 0.01
    sim = smart.interface.Smart(tag = 'highw_strato')
    sim.set_planet_proxima_b()    
    r_km = 149598000 * sim.smartin.r_AU

    
    fig, ax = plt.subplots(figsize = (30, 10))
    temp = avg_flux / avg_sflux 
    fpfs = avg_flux/avg_sflux * (sim.smartin.radius/r_km)**2
    
    ax.plot(wl[:len(flux)], temp[:len(wl)])
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$m)")
    ax2 = ax.twinx()
    ax2.set_ylabel("Planet-to-star contrast ratio")     
    ax2.plot(wl[:len(flux)], fpfs[:len(wl)])
    label = "Simulated Earth-like planet orbiting Proxima Centauri"
    ax.set_title(label)
    ax.set_xlim(0.5,2)
    ax.axvspan(0.61, 0.65, alpha=0.5, color='0.85')
    ax.axvspan(0.67, 0.71, alpha=0.5, color='0.85')
    ax.axvspan(0.74, 0.78, alpha=0.5, color='0.85')
    ax.axvspan(1.25, 1.29, alpha=0.5, color='0.85')
    fig.savefig("/gscratch/vsm/mwjl/projects/high_res/plots/prox_newCIA.png", bbox_inches = 'tight')
    
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
        final_plot()
    else:
        final_plot()
