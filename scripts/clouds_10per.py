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

def clouds(res, cirrus, strato, lamin, lamax):
    
    sim = smart.interface.Smart(tag = "earth_10per")
    sim.smartin.alb_file = "composite1_txt.txt"
    infile = "earth_avg_10per.pt"

    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, "clouds")

    try:
        os.mkdir(place)
    except OSError:
        pass

        
    sim.set_run_in_place(place) 
    sim.set_executables_automatically()

    sim.load_atmosphere_from_pt(infile, addn2 = False)
    o2 = sim.atmosphere.gases[6]
    o2.cia_file = "cia_adj_calc_10per.cia"
    o2.xsec_file = "dat_10per.dat"


    sim.smartin.FWHM = res
    sim.smartin.sample_res = res

    sim.smartin.minwn = 1e4/lamax
    sim.smartin.maxwn = 1e4/lamin 

    sim.lblin.minwn = 1e4/lamax
    sim.lblin.maxwn = 1e4/lamin 


    sim.gen_lblscripts()
    sim.run_lblabc()

    if cirrus == 1:
        sim.aerosols = smart.interface.Aerosols(cirrus=True, stratocum=False)
        sim.tag = "_cirrus"

    elif strato == 1:
        sim.aerosols = smart.interface.Aerosols(cirrus=False, stratocum=True)
        sim.tag = "_strato"

    else:
        pass

    import platform
    
    if platform.system() == 'Darwin':
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

    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux

    adj_flux = flux/sflux

    return(wl, adj_flux)

def longplot(lamin, lamax):
    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, "clouds")

    try:
        os.mkdir(place)
    except OSError:
        pass
    
    res = 0.01
    
    sim = smart.interface.Smart(tag = "prox")
    infile = "earth_avg_10per.pt"
    label = "Earth"
    sim.smartin.alb_file = "composite1_txt.txt"
    sim.set_run_in_place(place) 
    sim.set_executables_automatically()

    sim.load_atmosphere_from_pt(infile, addn2 = False)
    o2 = sim.atmosphere.gases[6]
    o2.cia_file = "cia_adj_calc_10per.cia"

    sim.smartin.sza = 57

    sim.smartin.FWHM = res
    sim.smartin.sample_res = res

    sim.smartin.minwn = 1e4/lamax
    sim.smartin.maxwn = 1e4/lamin 

    sim.lblin.minwn = 1e4/lamax
    sim.lblin.maxwn = 1e4/lamin 


    sim.gen_lblscripts()
    sim.run_lblabc()
    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux

    adj_flux = flux/sflux
    
    return(wl, adj_flux)


def plotting(lamin, lamax):
    cirrus_wl, cirrus_flux = clouds(0.01, 1, 0, lamin, lamax)
    strato_wl, strato_flux = clouds(0.01, 0, 1, lamin, lamax)
    wl, flux = longplot(lamin, lamax)
    length_wl = min(len(cirrus_wl), len(strato_wl), len(wl))-1
    avg_wl = (cirrus_wl[:length_wl] + strato_wl[:length_wl] + wl[:length_wl])/3
    avg_flux = (cirrus_flux[:length_wl] + strato_flux[:length_wl] +flux[:length_wl])/3

    fig_name = int(100*(float(lamin) + float(lamax))/2)

    fig, ax = plt.subplots(figsize = (30, 10))
    ax3 = ax.twinx()
    ax3.plot(avg_wl, avg_flux)
    ax3.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$ m)")
    ax.set_title("Earth Plot with 10% Oxygen")
    ax.plot(avg_wl, avg_flux)
    fig.savefig(str(fig_name) + "_10per.png", bbox_inches = 'tight')
    
if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="tenper",
                               subname="submit.csh",
                               workdir = "",
                               nodes = 1,
                               mem = "500G",
                               walltime = "48:00:00",
                               ntasks = 28,
                               account = "vsm",
                               submit = True,
                               rm_after_submit = True)
    elif platform.node().startswith("n"):
        # On a mox compute node: ready to run
        plotting(0.61,0.65)
        plotting(0.67,0.71)
        plotting(0.74,0.78)
        plotting(1.25,1.29)
    else:
        # Presumably, on a regular computer: ready to run
        plotting(0.61,0.65)
        plotting(0.67,0.71)
        plotting(0.74,0.78)
        plotting(1.25,1.29)







