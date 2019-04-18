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

def clouds(res, cirrus, strato):
    lamin = 0.5
    lamax = 2.0
    
    sim = smart.interface.Smart(tag = "earth")
    sim.smartin.alb_file = "composite1_txt.txt"
    infile = "earth_avg.pt"

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
    o2.cia_file = "cia_adj_calc.cia"
    o2.xsec_file = "dat_10per.dat"

    sim.smartin.FWHM = res
    sim.smartin.sample_res = res

    sim.smartin.minwn = 1e4/lamax
    sim.smartin.maxwn = 1e4/lamin 

    sim.lblin.minwn = 1e4/lamax
    sim.lblin.maxwn = 1e4/lamin 


    sim.gen_lblscripts()
    sim.run_lblabc()

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
        clouds(0.01, 1, 1):
    else:
        # Presumably, on a regular computer: ready to run
        clouds(0.01, 1, 1):








