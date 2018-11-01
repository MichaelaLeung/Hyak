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

def longplot(atmos, res, lamin, lamax, cirrus, strato, counter):
    sim = smart.interface.Smart(tag = str(atmos)+"sensi")
    sim.smartin.alb_file = "composite1_txt.txt"
    if atmos == "earth":
        infile = "earth_avg.pt"
    elif atmos == "prox":
        infile = "profile_Earth_proxb_.pt_filtered"
        label = "Self Consistent PCb"
        sim.set_planet_proxima_b()
    elif atmos == "highd":
        infile = "10bar_O2_dry.pt_filtered.pt"
        label = "10 bar O2 PCb"
        sim.set_planet_proxima_b()
    elif atmos == "highw":
        infile = "10bar_O2_wet.pt_filtered.pt"
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


    sim.load_atmosphere_from_pt(infile, addn2 = False)


    sim.smartin.FWHM = res
    sim.smartin.sample_res = res

    sim.smartin.minwn = 1e4/lamax
    sim.smartin.maxwn = 1e4/lamin 

    sim.lblin.minwn = 1e4/lamax
    sim.lblin.maxwn = 1e4/lamin
    sim.run_gas_sensitivity_test(run=True)

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

    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux

    adj_flux = flux/sflux * ((sim.smartin.radius / sim.smartin.r_AU) **2 )

if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="sensi",
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
        longplot("highd", 0.01, 0.5, 2, False, False, 0)
        counter = 0
    else:
        longplot("highd", 10, 0.6, 1.3, False, False, 0)
        counter = 0








