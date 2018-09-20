import numpy as np
import smart
from matplotlib import pyplot as plt
import matplotlib
from matplotlib.collections import LineCollection
from astropy.io import fits 
from matplotlib.collections import LineCollection
import matplotlib
matplotlib.rcParams['text.usetex'] = False
import sys, os

def integration_metric(lamin, lamax, mode):
    res = 0.01
    res2 = 1
    tag = "hipass"
    sim = smart.interface.Smart(tag = tag)
    if mode == 0: 
        infile = "profile_Earth_proxb_.pt_filtered"
        info = "reg"
    else:
        infile = "10bar_O2_dry.pt_filtered.pt"
        info = "high_o2"
    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, "integration")
    try:
        os.mkdir(place)
    except OSError:
            pass        

    sim.set_run_in_place(place) 
    sim.set_executables_automatically()
    sim.set_planet_proxima_b()
    sim.set_star_proxima()
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


    tag3 = tag + "_low"
    sim3 = smart.interface.Smart(tag = tag3)
    sim3.set_run_in_place(place)
    
    sim3.set_executables_automatically()
    sim.set_planet_proxima_b()
    sim.set_star_proxima()    sim3.load_atmosphere_from_pt(infile, addn2 = False)
    sim3.smartin.FWHM = res2
    sim3.smartin.sample_res = res2
    sim3.smartin.minwn = 1e4/lamax
    sim3.smartin.maxwn = 1e4/lamin 
    sim3.lblin.minwn = 1e4/lamax
    sim3.lblin.maxwn = 1e4/lamin 
    sim3.gen_lblscripts()
    sim3.run_lblabc()
    sim3.write_smart()
    sim3.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux

    flux = flux/sflux

    sim3.open_outputs()
    wl_low = sim3.output.rad.lam
    flux_low = sim3.output.rad.pflux
    sflux_low = sim3.output.rad.sflux

    flux_low = flux_low/sflux_low

    long_flux = []
    for i in flux_low:
        j = 0
        while j < 100: 
            long_flux.append(i)
            j = j+1

    mixed = []
    i = 0
    while i < len(flux):
        temp = (flux[i] + long_flux[i]) / 2
        mixed.append(temp)
        i = i+1


    i = 0
    flattened = []
    while i < len(flux)- 25: 
        avg = np.mean(flux[i:i+25])
        j = 0
        while j < 25:
            flattened.append(avg)
            j = j+1
        i = i+25
        

    out = []
    i = 0
    while i < len(mixed[:-25]): 
        diff = abs(mixed[i] - flattened[i])
        out.append(diff)
        i = i+1


    import scipy.integrate as integrate
    adds = integrate.trapz(out, wl[:-25])
    name = str(lamin) + "to" + str(lamax), str(info), str(abs(adds))
    f = open("integrations.txt", "a")
    f.write(str(name) + "\n")


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
        integration_metric(0.61,0.65,0)
        integration_metric(0.66,0.7,0)
        integration_metric(0.74, 0.78,0)
        integration_metric(0.61,0.65,1)
        integration_metric(0.66,0.7,1)
        integration_metric(0.74, 0.78,1)
    else:
        # Presumably, on a regular computer: ready to run
        integration_metric(0.62,0.64,0)

 
