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
import random

def ocean_loss(lamin, lamax):
    res = 1/(10*lamin)

    sim2 = smart.interface.Smart(tag = "highd")
    infile2 = "10bar_O2_dry.pt_filtered.pt"
    label = "Ocean Loss"
    sim2.smartin.alb_file = "desert_highd.alb"
    sim2.set_planet_proxima_b()
    sim2.set_star_proxima()

    sim2.set_run_in_place() 
    sim2.set_executables_automatically()

    sim2.smartin.sza = 57
    sim2.load_atmosphere_from_pt(infile2, addn2 = False, scaleP = 1.0)

    sim2.smartin.FWHM = res
    sim2.smartin.sample_res = res

    sim2.smartin.minwn = 1e4/lamax
    sim2.smartin.maxwn = 1e4/lamin 

    sim2.lblin.minwn = 1e4/lamax
    sim2.lblin.maxwn = 1e4/lamin 

    sim2.gen_lblscripts()
    sim2.run_lblabc()
    sim2.write_smart(write_file = True)
    sim2.run_smart()

    sim2.open_outputs()
    wl2 = sim2.output.rad.lam
    flux2 = sim2.output.rad.pflux
    sflux2 = sim2.output.rad.sflux

    adj_flux2 = flux2/sflux2
    return(wl2, adj_flux2)

def ocean_outgassing(lamin, lamax):
    res = 1/(10*lamin)

    sim2 = smart.interface.Smart(tag = "highw")
    infile2 = "10bar_O2_wet.pt_filtered.pt"
    label = "Ocean Outgassing"
    sim2.smartin.alb_file = "earth_noveg_highw.alb"
    sim2.set_planet_proxima_b()
    sim2.set_star_proxima()

    sim2.set_run_in_place() 
    sim2.set_executables_automatically()

    sim2.smartin.sza = 57
    sim2.load_atmosphere_from_pt(infile2, addn2 = False, scaleP = 1.0)

    sim2.smartin.FWHM = res
    sim2.smartin.sample_res = res

    sim2.smartin.minwn = 1e4/lamax
    sim2.smartin.maxwn = 1e4/lamin 

    sim2.lblin.minwn = 1e4/lamax
    sim2.lblin.maxwn = 1e4/lamin 

    sim2.gen_lblscripts()
    sim2.run_lblabc()
    sim2.write_smart(write_file = True)
    sim2.run_smart()

    sim2.open_outputs()
    wl2 = sim2.output.rad.lam
    flux2 = sim2.output.rad.pflux
    sflux2 = sim2.output.rad.sflux

    adj_flux2 = flux2/sflux2
    return(wl2, adj_flux2)

def ocean_loss_noO4(lamin, lamax):
    res = 1/(10*lamin)

    sim2 = smart.interface.Smart(tag = "highd")
    infile2 = "10bar_O2_dry.pt_filtered.pt"
    label = "Ocean Loss"
    sim2.smartin.alb_file = "desert_highd.alb"
    sim2.set_planet_proxima_b()
    sim2.set_star_proxima()

    sim2.set_run_in_place() 
    sim2.set_executables_automatically()

    sim2.smartin.sza = 57
    sim2.load_atmosphere_from_pt(infile2, addn2 = False, scaleP = 1.0)

    sim2.smartin.FWHM = res
    sim2.smartin.sample_res = res

    sim2.smartin.minwn = 1e4/lamax
    sim2.smartin.maxwn = 1e4/lamin 

    sim2.lblin.minwn = 1e4/lamax
    sim2.lblin.maxwn = 1e4/lamin

    o2 = sim2.atmosphere.gases[1]
    o2.cia_file = None

    sim2.gen_lblscripts()
    sim2.run_lblabc()
    sim2.write_smart(write_file = True)
    sim2.run_smart()

    sim2.open_outputs()
    wl2 = sim2.output.rad.lam
    flux2 = sim2.output.rad.pflux
    sflux2 = sim2.output.rad.sflux

    adj_flux2 = flux2/sflux2
    return(wl2, adj_flux2)

def ocean_outgassing_noO4(lamin, lamax):
    res = 1/(10*lamin)

    sim2 = smart.interface.Smart(tag = "highw")
    infile2 = "10bar_O2_wet.pt_filtered.pt"
    label = "Ocean Outgassing"
    sim2.smartin.alb_file = "earth_noveg_highw.alb"
    sim2.set_planet_proxima_b()
    sim2.set_star_proxima()

    sim2.set_run_in_place() 
    sim2.set_executables_automatically()

    sim2.smartin.sza = 57
    sim2.load_atmosphere_from_pt(infile2, addn2 = False, scaleP = 1.0)

    sim2.smartin.FWHM = res
    sim2.smartin.sample_res = res

    sim2.smartin.minwn = 1e4/lamax
    sim2.smartin.maxwn = 1e4/lamin 

    sim2.lblin.minwn = 1e4/lamax
    sim2.lblin.maxwn = 1e4/lamin

    o2 = sim2.atmosphere.gases[1]
    o2.cia_file = None

    sim2.gen_lblscripts()
    sim2.run_lblabc()
    sim2.write_smart(write_file = True)
    sim2.run_smart()

    sim2.open_outputs()
    wl2 = sim2.output.rad.lam
    flux2 = sim2.output.rad.pflux
    sflux2 = sim2.output.rad.sflux

    adj_flux2 = flux2/sflux2
    return(wl2, adj_flux2)



def plotting(lamin, lamax, atmos, title):
    fig_name = int(100*(float(lamin) + float(lamax))/2)
    if atmos == 0: #0 = ocean loss
        wl, flux = ocean_loss(lamin, lamax)
        wl2, flux2 = ocean_loss_noO4(lamin, lamax)
        fig, ax = plt.subplots(figsize = (10,10))
        ax.plot(wl, flux, label = "Ocean Loss")
        ax.plot(wl2, flux2, label = "Ocean Loss, no O2O2")
        ax.set_title(title)
        ax.set_ylabel("Reflectance")
        ax.set_xlabel("Wavelength ($\mu$ m)")
        ax.legend()
        fig.savefig(str(fig_name) +  "noO4.png", bbox_inches = "tight")
    else:
        wl, flux = ocean_outgassing(lamin, lamax)
        wl2, flux2 = ocean_outgassing_noO4(lamin, lamax)
        fig, ax = plt.subplots(figsize = (10,10))
        ax.plot(wl, flux, label = "Ocean Outgassing")
        ax.plot(wl2, flux2, label = "Ocean Outgassing, no O2O2")
        ax.set_title(title)
        ax.set_ylabel("Reflectance")
        ax.set_xlabel("Wavelength ($\mu$ m)")
        ax.legend()
        fig.savefig(str(fig_name) +  "noO4_ocean.png", bbox_inches = "tight")

   
if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="smartplt",
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
        plotting(0.61,0.65,0,"Gamma band (0.63) Ocean Loss")
        plotting(0.67,0.71,0,"Oxygen B band (0.69) Ocean Loss")
        plotting(0.74,0.78,0,"Oxygen A band (0.76) Ocean Loss")
        plotting(1.25,1.29,0,"1.27 Ocean Loss")
        plotting(0.61,0.65,1,"Gamma band (0.63) Ocean Outgassing")
        plotting(0.67,0.71,1,"Oxygen B band (0.69) Ocean Outgassing")
        plotting(0.74,0.78,1,"Oxygen A band (0.76) Ocean Outgassing")
        plotting(1.25,1.29,1,"1.27 Ocean Loss")
    else:
        plotting(0.74,0.78,0,"0.76 Ocean Loss")
