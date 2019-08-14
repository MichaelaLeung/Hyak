import numpy as np
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from astropy.io import fits
import smart
import sys, os
import datetime
matplotlib.rcParams['text.usetex'] = False
import random

def clouds(res, cirrus, strato):
    lamin = 0.5
    lamax = 2.0
    
    sim = smart.interface.Smart(tag = "prox")
    sim.smartin.alb_file = "composite1_txt.txt"
    infile = "profile_Earth_proxb_.pt_filtered"

    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, "clouds")

    try:
        os.mkdir(place)
    except OSError:
        pass

        
    sim.set_run_in_place(place) 
    sim.set_executables_automatically()
    
    sim.lblin.par_file = '/gscratch/vsm/alinc/fixed_input/HITRAN2016' #/gscratch/vsm/alinc/fixed_input/
    sim.lblin.hitran_tag = 'hitran2016'
    sim.lblin.fundamntl_file = '/gscratch/vsm/alinc/fixed_input/fundamntl2016.dat'
    sim.lblin.lblabc_exe = '/gscratch/vsm/alinc/exec/lblabc_2016'

    sim.load_atmosphere_from_pt(infile, addn2 = False)
    o2 = sim.atmosphere.gases[6]
    o2.cia_file = "cia_adj_mix.cia"

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
    
    if platform.system() == 'Jarvis':
        # On a Mac: usetex ok
        matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})
        matplotlib.rcParams['font.size'] = 25.0
        matplotlib.rc('text', usetex=True)
    elif platform.node().startswith("n"):
        # On hyak: usetex not ok, must change backend to 'agg'
        matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})
        matplotlib.rcParams['font.size'] = 25.0
        matplotlib.rc('text', usetex=False)
        plt.switch_backend('agg')

    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux

    adj_flux = flux/sflux

    return(wl, adj_flux)

def longplot():
    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, "clouds")

    try:
        os.mkdir(place)
    except OSError:
        pass
    
    lamin = 0.5
    lamax = 2.0 
    
    res = 0.01
    
    sim = smart.interface.Smart(tag = "prox")
    infile = "profile_Earth_proxb_.pt_filtered"
    label = "Prox"
    sim.smartin.alb_file = "composite1_txt.txt"
    sim.set_run_in_place(place) 
    sim.set_executables_automatically()

    sim.lblin.par_file = '/gscratch/vsm/alinc/fixed_input/HITRAN2016' #/gscratch/vsm/alinc/fixed_input/
    sim.lblin.hitran_tag = 'hitran2016'
    sim.lblin.fundamntl_file = '/gscratch/vsm/alinc/fixed_input/fundamntl2016.dat'
    sim.lblin.lblabc_exe = '/gscratch/vsm/alinc/exec/lblabc_2016'

    sim.load_atmosphere_from_pt(infile, addn2 = False)
    o2 = sim.atmosphere.gases[6]
    o2.cia_file = "cia_adj_mix.cia"

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


def plotting():
    cirrus_wl, cirrus_flux = clouds(0.01, 1, 0)
    strato_wl, strato_flux = clouds(0.01, 0, 1)
    wl, flux = longplot()
    length_wl = min(len(cirrus_wl), len(strato_wl), len(wl))-1
    avg_wl = (cirrus_wl[:length_wl] + strato_wl[:length_wl] + wl[:length_wl])/3
    avg_flux = (cirrus_flux[:length_wl] + strato_flux[:length_wl] +flux[:length_wl])/3
    fig, ax = plt.subplots(figsize = (30, 10))
    ax3 = ax.twinx()
    ax3.plot(avg_wl, avg_flux)
    ax3.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$ m)")
    ax.set_title("Simulated Earth-like planet with clouds orbiting Proxima Centauri")
    ax.set_xlim(0.5,2)
    ax.axvspan(0.61, 0.65, alpha=0.5, color='0.85')
    ax.axvspan(0.67, 0.71, alpha=0.5, color='0.85')
    ax.axvspan(0.74, 0.78, alpha=0.5, color='0.85')
    ax.axvspan(1.25, 1.29, alpha=0.5, color='0.85')
    ax.plot(avg_wl, avg_flux)
    fig.savefig("avg_clouds.png", bbox_inches = 'tight')
    
if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="clouds",
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
        plotting()
    else:
        # Presumably, on a regular computer: ready to run
        plotting()        








