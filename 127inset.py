

import numpy as np
import smart
from matplotlib import pyplot as plt
import matplotlib
from matplotlib.collections import LineCollection
from astropy.io import fits 
import coronagraph as cg
from pylab import *
import matplotlib
matplotlib.rcParams['text.usetex'] = False
import os

def run_smart(res, infile, lamin, lamax):
        sim = smart.interface.Smart(tag = "127inset")

        tag = "proxima_higho2_d"
        HERE = os.path.dirname(os.path.abspath(__file__))
        place = os.path.join(HERE, "127inset")
        try:
                os.mkdir(place)
        except OSError:
                pass        

        sim.set_run_in_place(place) 
        sim.set_executables_automatically()
 
        sim.set_planet_proxima_b()


        sim.load_atmosphere_from_pt(infile, addn2 = False)

        sim.smartin.FWHM = res
        sim.smartin.sample_res = res
        sim.smartin.minwn = 1e4/lamax
        sim.smartin.maxwn = 1e4/lamin 

        sim.lblin.minwn = 1e4/lamax
        sim.lblin.maxwn = 1e4/lamin 

        sim.gen_lblscripts()
        sim.run_lblabc()

        sim.write_smart()

        sim.run_smart()

        sim.open_outputs()
        wl = sim.output.rad.lam
        flux = sim.output.rad.pflux
        return(wl, flux)


def inset(infile):
        wl, flux = run_smart(0.01, infile, 0.5, 2.0)
        wl_small, flux_small = run_smart(0.0787, infile, 1.24, 1.30)
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

        x = [1.24,1.24,1.30,1.30]

        y = [0, max(flux), max(flux), 0]

        fig, ax = plt.subplots(2, 1, figsize = (15,8))
        ax[0].plot(wl,flux)
        ax[0].fill(x,y, '0.75')
        ax[1].plot(wl_small, flux_small)
        name = infile
        ax[0].set_title(name)
        fig.savefig(name + ".png", bbox_inches = 'tight')




if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="127inset",
                               subname="submit.csh",
                               workdir = "",
                               nodes = 1,
                               mem = "500G",
                               walltime = "0",
                               ntasks = 28,
                               account = "vsm",
                               submit = True,
                               rm_after_submit = True)
    elif platform.node().startswith("n"):
        # On a mox compute node: ready to run
        inset('10bar_O2_dry.pt_filtered.pt')
        inset("profile_Earth_proxb_.pt_filtered")

    else:
        # Presumably, on a regular computer: ready to run
        infile = '10bar_O2_dry.pt_filtered.pt'
        wl, flux = run_smart(1, infile, 1.25, 1.32)
        wl_small, flux_small = run_smart(1, infile, 1.27, 1.29)

        x = [1.24,1.24,1.30,1.30]

        y = [0, max(flux), max(flux), 0]

        fig, ax = plt.subplots(2, 1, figsize = (15,8))
        ax[0].plot(wl,flux)
        ax[0].fill(x,y, '0.75')
        ax[1].plot(wl_small, flux_small)
        name = infile
        ax[0].set_title(name)
        fig.savefig(name + ".png", bbox_inches = 'tight')
        



