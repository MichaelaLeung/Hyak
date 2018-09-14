import numpy as np
import smart
from matplotlib import pyplot as plt
import matplotlib
from matplotlib.collections import LineCollection
from astropy.io import fits 
import coronagraph as cg
import matplotlib
matplotlib.rcParams['text.usetex'] = False

def longplot(atmos, res, lamin, lamax, cirrus, strato):
    
    sim = smart.interface.Smart(tag = atmos)
    if atmos = "earth":
        infile = "earth_avg.pt"
    elif atmos = "prox":
        infile = "profile_Earth_proxb_.pt_filtered"
    elif atmos = "10bard":
        infile = "10bar_O2_dry.pt_filtered.pt"
    elif atmos = "10barw":
        infile = "10bar_O2_wet.pt_filtered.pt"
    elif atmos = "arch_prox":
        infile = "clearsky_archean.pt"
        
    sim.set_run_in_place("/Users/mwl/python/") ** REPLACE THIS**
    sim.set_run_on_mac() ** REPLACE THIS** 
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

    if cirrus:
        sim.aerosols = smart.interface.Aerosols(cirrus=True, stratocum=False)
        sim.tag = atmos + "_cirrus"

    if strato:
        sim.aerosols = smart.interface.Aerosols(cirrus=False, stratocum=True)
        sim.tag = atmos + "_strato"


    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.read.sflux

    adj_flux = flux/sflux * ((sim.smartin.radius / sim.smartin.r_AU) **2 )

    fig, ax = plt.subplots(figsize = (30, 10))
    ax.plot(wl, flux)
    ax.set_title("earth: stratocumulus clouds")
    fig.savefig("earth_clouds_strat", bbox_inches = 'tight')

if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="longplt",
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
        longplot(earth, 0.01, 0.5, 2, True, False)
        longplot(earth, 0.01, 0.5, 2, False, True)
        longplot(prox, 0.01, 0.5, 2, False, False)
        longplot(arch_prox, 0.01, 0.5, 2, False, False)
        longplot(10bard, 0.01, 0.5, 2, False, False)
        longplot(10barw, 0.01, 0.5, 2, False, False)
    else:
        # Presumably, on a regular computer: ready to run
        longplot(earth, 1, 0.5, 0.52, True, False)
        longplot(earth, 1, 0.5, 0.52, False, True)
        longplot(prox, 1, 0.5, 0.52, False, False)
        longplot(arch_prox, 1, 0.5, 0.52, False, False)
        longplot(10bard, 1, 0.5, 0.52, False, False)
        longplot(10barw, 1, 0.5, 0.52, False, False)








