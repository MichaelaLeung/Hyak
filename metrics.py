import numpy as np
import smart
from matplotlib import pyplot as plt
import matplotlib
from matplotlib.collections import LineCollection
from astropy.io import fits 
import matplotlib
import sys, os
import datetime
matplotlib.rcParams['text.usetex'] = False
import scipy.integrate as integrate

def run_smart(lamin, lamax):
    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, "metrics")

    try:
        os.mkdir(place)
    except OSError:
        pass

        

    sim = smart.interface.Smart(tag = "prox")
    infile = "profile_Earth_proxb_.pt_filtered"
    res = 0.01
    low_res = 1
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

    sim2 = smart.interface.Smart(tag = "prox")

    sim2.set_run_in_place(place) 
    sim2.set_executables_automatically()
    sim2.load_atmosphere_from_pt(infile, addn2 = False)
    sim2.smartin.FWHM = low_res
    sim2.smartin.sample_res = low_res
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
    sflux = sim.output.rad.sflux
    adj_flux = flux/sflux * ((sim.smartin.radius / (sim.smartin.r_AU * 149598000)) **2 )

    sim2.open_outputs()
    wl_low = sim2.output.rad.lam
    flux_low = sim2.output.rad.pflux

    return (wl, flux, adj_flux, wl_low, flux_low)
 
    
def interval(wl, flux):
    lis = []
    st_wl = wl[0]
    delta = []
    exo_dict = dict(zip(flux, wl))
    final_dict = dict()
    cutoff = np.median(flux)
    for i in flux: 
        temp = (i > cutoff)
        if temp == True: 
            start = i 
            st_wl = exo_dict[start]
        elif temp == False: 
            finish = i 
            fin_wl = exo_dict[finish]
            delta = (fin_wl - st_wl)
            final_dict[st_wl] = fin_wl
    return(len(final_dict))

def high_pass(flux, flux_low):
    k = 0
    x = 0
    z = 0
    mixed = []
    out = []
    flattened = []
    long_flux = []
    for i in flux_low:
        j = 0
        while j < 100: 
            long_flux.append(i)
            j = j+1
    print(len(long_flux))
        
    while k < len(flux):
        temp = (flux[k] + long_flux[k]) / 2
        mixed.append(temp)
        k = k+1
    
    print(len(mixed))
        
    while x < len(flux)- 25: 
        avg = np.mean(flux[x:x+25])
        y = 0
        while y < 25:
            flattened.append(avg)
            y = y+1
        x = x+25
        
        
    while z < len(mixed[:-25]): 
        diff = abs(mixed[z] - flattened[z])
        out.append(diff)
        z = z+1
        
    return(out)


def fourier(flux):
    from scipy.fftpack import fft, rfft, fftfreq
    yf = rfft(flux)
    return(yf)


def outputs(lamin, lamax):
    wl, flux, adj_flux, wl_low, flux_low = run_smart(lamin, lamax)
    adds = max(abs(integrate.trapz((high_pass(flux, flux_low), wl[:-25]))))
    high = interval(wl[:-25], (high_pass(flux, flux_low)))
    fouri = interval(wl,fourier(flux))
    label = str(lamin) + "to" + str(lamax)
    out = label, "fpfs", np.median(adj_flux), "line cutoff", high, "integral", adds, "together", (np.median(adj_flux)*high*adds)
    f = open("outputs_small.txt", "a")
    f.write(str(out) + "\n")
    

if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="metrics",
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
        number = range(50,250, 2)
        for i in number:
            i = float(i)
            i = i/100
            outputs(i, i+0.1)
    else:
        # Presumably, on a regular computer: ready to run
       number = range(50,250, 2)
       for i in number:
            i = float(i)
            i = i/100
            outputs(i, i+0.1)
