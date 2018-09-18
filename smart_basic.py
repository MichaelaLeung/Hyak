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

def smart_basic(res, lamin, lamax):

    sim = smart.interface.Smart(tag = "prox")
    sim2 = smart.interface.Smart(tag = "earth")
    infile = "profile_Earth_proxb_.pt_filtered"
    earth_infile = "earth_avg.pt"
    
    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, "smart_basic")

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
    sim.gen_lblscripts()
    sim.run_lblabc()
    sim.write_smart(write_file = True)
    sim.run_smart()

    sim2.set_run_in_place(place) 
    sim2.set_executables_automatically()
    sim2.load_atmosphere_from_pt(infile, addn2 = False)
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

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux

    sim2.open_outputs()
    earth_wl = sim2.output.rad.lam
    earth_flux = sim2.output.rad.pflux


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

  #  fig, ax = plt.subplots(figsize=(12,10))
 #   plt.plot(wl,flux)
    
 

if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="manyplt",
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
        number = range(80,82, 1)
        for i in number:
            i = float(i)
            i = i/100
            smart_basic(1, i, i+0.1)
    else:
        # Presumably, on a regular computer: ready to run
        number = range(80,82,1)
        for i in number:
            i = float(i)
            i = i/100
            smart_basic(1, i, i+0.1)









