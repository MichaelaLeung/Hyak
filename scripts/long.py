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

def longplot(atmos, res, lamin, lamax, cirrus, strato):
    
    sim = smart.interface.Smart(tag = atmos)
    if atmos == "earth":
        infile = "earth_avg.pt"
        sim.smartin.alb_file = "composite1_txt.txt"
    elif atmos == "prox":
        infile = "profile_Earth_proxb_.pt_filtered"
        label = "Self Consistent PCb (Earth like)"
        sim.smartin.alb_file = "composite1_txt.txt"
        sim.set_planet_proxima_b()
    elif atmos == "highd":
        infile = "10bar_O2_dry.pt_filtered.pt"
        label = "10 bar O2 PCb"
        sim.smartin.alb_file = "desert_highd.alb"
        sim.set_planet_proxima_b()
    elif atmos == "highw":
        infile = "10bar_O2_wet.pt_filtered.pt"
        label = "10 bar O2 PCb with water vapor"
        sim.smartin.alb_file = "earth_noveg_highw.alb"
        sim.set_planet_proxima_b()
    elif atmos == "arch_prox":
        infile = "clearsky_archean.pt"
        sim.set_planet_proxima_b()
    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, "longplot")

    try:
        os.mkdir(place)
    except OSError:
        pass

        
    sim.set_run_in_place(place) 
    sim.set_executables_automatically()

    sim.smartin.sza = 57
    sim.load_atmosphere_from_pt(infile, addn2 = False)


    sim.smartin.FWHM = res
    sim.smartin.sample_res = res

    sim.smartin.minwn = 1e4/lamax
    sim.smartin.maxwn = 1e4/lamin 

    sim.lblin.minwn = 1e4/lamax
    sim.lblin.maxwn = 1e4/lamin 


    sim.gen_lblscripts()
    sim.run_lblabc()

    if cirrus == True:
        sim.aerosols = smart.interface.Aerosols(cirrus=True, stratocum=False)
        sim.tag = atmos + "_cirrus"

    elif strato == True:
        sim.aerosols = smart.interface.Aerosols(cirrus=False, stratocum=True)
        sim.tag = atmos + "_strato"

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

    fig, ax = plt.subplots(figsize = (30, 10))
    ax.plot(wl, adj_flux)
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$ m)")
    ax.set_title(label)
    fig.savefig(str(atmos) + ".png", bbox_inches = 'tight')




if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="long",
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
        longplot("earth", 0.01, 0.5, 2, True, False)
 #       longplot("earth", 0.01, 0.5, 2, False, True)
 #       longplot("prox", 0.01, 0.5, 2, False, False)
 #       longplot("highd", 0.01, 0.5, 2, False, False)
 #       longplot("highw", 0.01, 0.5, 2, False, False)
 #       longplot("arch_prox", 0.01, 0.5, 2, False, False)
    else:
        # Presumably, on a regular computer: ready to run
        longplot("earth", 1, 0.5, 0.501, True, False)
        longplot("earth", 1, 0.5, 0.501, False, True)
        longplot("prox", 1, 0.5, 0.501, False, False)
        longplot("highd", 1, 0.5, 0.501, False, False)
        longplot("highw", 1, 0.5, 0.501, False, False)
        longplot("arch_prox", 1, 0.5, 0.501, False, False)








