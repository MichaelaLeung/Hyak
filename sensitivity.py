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
        label = "Self Consistent PCb (Earth like)"
        sim.set_planet_proxima_b()
    elif atmos == "highd":
        infile = "10bar_O2_dry.pt_filtered.pt"
        label = "10 bar O2 PCb"
        sim.set_planet_proxima_b()
    elif atmos == "highw":
        infile = "10bar_O2_wet.pt_filtered.pt"
        label = "10 bar O2 PCb with water vapor"
        sim.set_planet_proxima_b()
    elif atmos == "arch_prox":
        infile = "clearsky_archean.pt"
        sim.set_planet_proxima_b()
    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, "sensi")

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
    sim.run_gas_sensitivity_test(run=True)

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
        longplot("highw", 0.01, 0.5,0.6, False, False, 0)
        longplot("highw", 0.01, 0.6,0.7, False, False, 0)
        longplot("highw", 0.01, 0.7,0.8, False, False, 0)
        longplot("highw", 0.01, 0.8,0.9, False, False, 0)
        longplot("highw", 0.01, 0.9,1.0, False, False, 0)
        longplot("highw", 0.01, 1.0,1.1, False, False, 0)
        longplot("highw", 0.01, 1.1,1.2, False, False, 0)
        longplot("highw", 0.01, 1.2,1.3, False, False, 0)
        longplot("highw", 0.01, 1.3,1.4, False, False, 0)
        longplot("highw", 0.01, 1.4,1.5, False, False, 0)
        longplot("highw", 0.01, 1.5,1.6, False, False, 0)
        longplot("highw", 0.01, 1.6,1.7, False, False, 0)
        longplot("highw", 0.01, 1.7,1.8, False, False, 0)
        longplot("highw", 0.01, 1.8,1.9, False, False, 0)
        longplot("highw", 0.01, 1.9,2.0, False, False, 0)
        counter = 0
    else:
        longplot("highw", 0.01, 0.8,0.9, False, False, 0)
        longplot("highw", 0.01, 0.9,1.0, False, False, 0)
        longplot("highw", 0.01, 1.0,1.1, False, False, 0)
        longplot("highw", 0.01, 1.1,1.2, False, False, 0)
        longplot("highw", 0.01, 1.2,1.3, False, False, 0)
        longplot("highw", 0.01, 1.3,1.4, False, False, 0)

        counter = 0








