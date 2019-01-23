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

def smart_basic(lamin, lamax, title, atmos):
        
    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, "oxygen_outplot")

    try:
        os.mkdir(place)
    except OSError:
        pass
    
    if atmos == 'dry':
        infile1 = "10bar_O2_dry.pt_filtered.pt"
        sim1 = smart.interface.Smart(tag = info2)
        sim1.load_atmosphere_from_pt(infile1, addn2 = False)
        sim1.smartin.alb_file = "desert_highd.alb"
        o2 = sim1.atmosphere.gases[1]
        o2.cia_file = None
        infile2 = "10bar_O2_dry.pt_filtered.pt"
        info2 = "highd"
        sim2 = smart.interface.Smart(tag = info2)
        sim2.load_atmosphere_from_pt(infile2, addn2 = False, scaleP = 1.0)
        sim2.smartin.alb_file = "desert_highd.alb"
    else:
        infile1 = "10bar_O2_wet.pt_filtered.pt"
        sim1 = smart.interface.Smart(tag = info2)
        sim1.load_atmosphere_from_pt(infile1, addn2 = False)
        sim1.smartin.alb_file = "earth_noveg_highw.alb"
        o2 = sim1.atmosphere.gases[1]
        o2.cia_file = None
        infile2 = "10bar_O2_wet.pt_filtered.pt"
        info2 = "highw"
        sim2 = smart.interface.Smart(tag = info2)
        sim2.load_atmosphere_from_pt(infile2, addn2 = False, scaleP = 1.0)
        sim2.smartin.alb_file = "earth_noveg_highw.alb"

    res = 1/(10*lamin)


    for sim in (sim1, sim2):
        sim.set_run_in_place()    
        sim.set_executables_automatically()
        sim.set_planet_proxima_b()
        sim.set_star_proxima()
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


    sim1.open_outputs()
    wl = sim1.output.rad.lam
    flux = sim1.output.rad.pflux
    sflux = sim1.output.rad.sflux
    flux = flux/sflux

    sim2.open_outputs()
    wl2 = sim2.output.rad.lam
    flux2 = sim2.output.rad.pflux
    sflux2 = sim2.output.rad.sflux
    flux2 = flux2/sflux2


    import platform
    if platform.system() == 'Jarvis':
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


    fig, ax = plt.subplots(figsize = (10,10))
    ax.plot(wl, flux, label = "Ocean Loss without O2-O2")
    ax.plot(wl2, flux2, label = "Ocean Loss")


    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$ m)")
    ax.set_title(title)
    fig_name = int(100*(float(lamin) + float(lamax))/2)
    ax.legend()
    if lamax > 1.29:
        ax.set_xlim(1.25,1.29)
    if atmos == 'dry':
        fig.savefig(str(fig_name) +  "_noO4.png", bbox_inches = "tight")
    else:
        fig.savefig(str(fig_name) +  "_noO4_ocean.png", bbox_inches = "tight")

 

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
        smart_basic(0.61, 0.65, "Gamma band O2-O2 CIA", 'dry')
        smart_basic(0.67, 0.71, "Oxygen B band O2-O2 CIA", 'dry')
        smart_basic(0.74, 0.78, "Oxygen A band O2-O2 CIA", 'dry')
        smart_basic(1.25,1.29, "1.27 band O2-O2 CIA", 'dry')
        smart_basic(0.61, 0.65, "Gamma band O2-O2 CIA", 'wet')
        smart_basic(0.67, 0.71, "Oxygen B band O2-O2 CIA", 'wet')
        smart_basic(0.74, 0.78, "Oxygen A band O2-O2 CIA", 'wet')
        smart_basic(1.25,1.29, "1.27 band O2-O2 CIA", 'wet')

    else:
        # Presumably, on a regular computer: ready to run
        smart_basic(0.61, 0.65, "Gamma band O2-O2 CIA", 'dry')
        smart_basic(0.67, 0.71, "Oxygen B band O2-O2 CIA", 'dry')
        smart_basic(0.74, 0.78, "Oxygen A band O2-O2 CIA", 'dry')
        smart_basic(1.25,1.29, "1.27 band O2-O2 CIA", 'dry')












