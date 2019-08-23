import numpy as np
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from astropy.io import fits
import smart
import sys, os
import datetime
matplotlib.rcParams['text.usetex'] = False
import random
import math

def clouds(lamin, lamax, cloud_type, res):

    place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim = smart.interface.Smart(tag = "prox")
    sim.set_run_in_place(place)
    
    sim.smartin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.lblin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.smartin.abs_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'

    infile = "/gscratch/vsm/mwjl/projects/high_res/inputs/profile_Earth_proxb_.pt_filtered"
    label = "Simulated Earth-like planet orbiting Proxima Centauri"
    sim.smartin.alb_file = "/gscratch/vsm/mwjl/projects/high_res/inputs/composite1_txt.txt"
    sim.set_planet_proxima_b()
    sim.load_atmosphere_from_pt(infile, addn2 = False)
    
    o2 = sim.atmosphere.gases[3]
    o2.cia_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/o4_calc.cia'
    label = "Earth-Like"
    sim.set_planet_proxima_b()
    sim.set_star_proxima()

    sim.set_executables_automatically()

    sim.lblin.par_file = '/gscratch/vsm/alinc/fixed_input/HITRAN2016.par' #/gscratch/vsm/alinc/fixed_input/
    sim.lblin.hitran_tag = 'hitran2016'
    sim.lblin.fundamntl_file = '/gscratch/vsm/alinc/fixed_input/fundamntl2016.dat'
    sim.lblin.lblabc_exe = '/gscratch/vsm/alinc/exec/lblabc_2016'
    sim.lblin.par_index = 7


    sim.smartin.sza = 57

    sim.smartin.FWHM = res
    sim.smartin.sample_res = res

    sim.smartin.minwn = 1e4/lamax
    sim.smartin.maxwn = 1e4/lamin 

    sim.lblin.minwn = 1e4/lamax
    sim.lblin.maxwn = 1e4/lamin 


    sim.gen_lblscripts()
    sim.run_lblabc()

    if cloud_type == 0:
        sim.aerosols = smart.interface.Aerosols(cirrus=True, stratocum=False)
        sim.tag = "prox_cirrus"

    elif cloud_type == 1:
        sim.aerosols = smart.interface.Aerosols(cirrus=False, stratocum=True)
        sim.tag = "prox_strato"

    else:
        pass

    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux

    adj_flux = flux/sflux * math.pi

    return(wl, adj_flux)

def earth_like(lamin, lamax, res, cia):
    place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim = smart.interface.Smart(tag = "prox")
    sim.set_run_in_place(place)
    infile = "/gscratch/vsm/mwjl/projects/high_res/inputs/profile_Earth_proxb_.pt_filtered"
    label = "Earth-Like"
    sim.smartin.alb_file = "/gscratch/vsm/mwjl/projects/high_res/inputs/composite1_txt.txt"
    sim.set_planet_proxima_b()
    sim.set_star_proxima()

    sim.smartin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.lblin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.smartin.abs_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'

    sim.set_executables_automatically()
    
    sim.lblin.par_file = '/gscratch/vsm/alinc/fixed_input/HITRAN2016.par' #/gscratch/vsm/alinc/fixed_input/
    sim.lblin.hitran_tag = 'hitran2016'
    sim.lblin.fundamntl_file = '/gscratch/vsm/alinc/fixed_input/fundamntl2016.dat'
    sim.lblin.lblabc_exe = '/gscratch/vsm/alinc/exec/lblabc_2016'
    sim.lblin.par_index	= 7

    sim.smartin.sza = 57
    sim.load_atmosphere_from_pt(infile, addn2 = False)

    sim.smartin.FWHM = res
    sim.smartin.sample_res = res

    sim.smartin.minwn = 1e4/lamax
    sim.smartin.maxwn = 1e4/lamin 

    sim.lblin.minwn = 1e4/lamax
    sim.lblin.maxwn = 1e4/lamin
    
    if cia == "new":
        o2 = sim.atmosphere.gases[3]
        o2.cia_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/o4_calc.cia'
    elif cia == "none":
        o2 = sim.atmosphere.gases[3]
        o2.cia_file = None
    else:
        pass 

    sim.gen_lblscripts()
    sim.run_lblabc()
    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux
    adj_flux = flux/sflux * math.pi

    cirrus_wl, cirrus_flux = clouds(lamin, lamax, 0, res)
    strato_wl, strato_flux = clouds(lamin, lamax,1, res)

    avg_flux = (flux + cirrus_flux + strato_flux)/3
    return(wl, avg_flux)

def ocean_loss(lamin, lamax, res, cia):
    place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim2 = smart.interface.Smart(tag = "highd")
    sim2.set_run_in_place(place) 
    infile2 = "/gscratch/vsm/mwjl/projects/high_res/inputs/10bar_O2_dry.pt_filtered.pt"
    label = "Ocean Loss"
    sim2.smartin.alb_file = "/gscratch/vsm/mwjl/projects/high_res/inputs/desert_highd.alb"
    sim2.set_planet_proxima_b()
    sim2.set_star_proxima()

    sim2.smartin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim2.lblin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim2.smartin.abs_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'

    sim2.set_executables_automatically()

    sim2.lblin.par_file = '/gscratch/vsm/alinc/fixed_input/HITRAN2016.par' #/gscratch/vsm/alinc/fixed_input/
    sim2.lblin.hitran_tag = 'hitran2016'
    sim2.lblin.fundamntl_file = '/gscratch/vsm/alinc/fixed_input/fundamntl2016.dat'
    sim2.lblin.lblabc_exe = '/gscratch/vsm/alinc/exec/lblabc_2016'
    sim2.lblin.par_index = 7

    sim2.smartin.sza = 57
    sim2.load_atmosphere_from_pt(infile2, addn2 = False, scaleP = 1.0)

    sim2.smartin.FWHM = res
    sim2.smartin.sample_res = res

    sim2.smartin.minwn = 1e4/lamax
    sim2.smartin.maxwn = 1e4/lamin 

    sim2.lblin.minwn = 1e4/lamax
    sim2.lblin.maxwn = 1e4/lamin

    if cia == "new":
        o2 = sim2.atmosphere.gases[1]
        o2.cia_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/o4_calc.cia'
    elif cia == "none":
        o2 = sim2.atmosphere.gases[1]
        o2.cia_file = None
    else:
        pass 

    sim2.gen_lblscripts()
    sim2.run_lblabc()
    sim2.write_smart(write_file = True)
    sim2.run_smart()

    sim2.open_outputs()
    wl2 = sim2.output.rad.lam
    flux2 = sim2.output.rad.pflux
    sflux2 = sim2.output.rad.sflux

    adj_flux2 = flux2/sflux2 * math.pi
    return(wl2, adj_flux2)

def ocean_outgassing(lamin, lamax, res, cia):
    place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim2 = smart.interface.Smart(tag = "highw")
    sim2.set_run_in_place(place) 
    infile2 = "/gscratch/vsm/mwjl/projects/high_res/inputs/10bar_O2_wet.pt_filtered.pt"
    label = "Ocean Outgassing"
    sim2.smartin.alb_file = "/gscratch/vsm/mwjl/projects/high_res/inputs/earth_noveg_highw.alb"
    sim2.set_planet_proxima_b()
    sim2.set_star_proxima()

    sim2.smartin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim2.lblin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim2.smartin.abs_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'

    sim2.set_executables_automatically()

    sim2.lblin.par_file = '/gscratch/vsm/alinc/fixed_input/HITRAN2016.par' #/gscratch/vsm/alinc/fixed_input/
    sim2.lblin.hitran_tag = 'hitran2016'
    sim2.lblin.fundamntl_file = '/gscratch/vsm/alinc/fixed_input/fundamntl2016.dat'
    sim2.lblin.lblabc_exe = '/gscratch/vsm/alinc/exec/lblabc_2016'
    sim2.lblin.par_index = 7

    sim2.smartin.sza = 57
    sim2.load_atmosphere_from_pt(infile2, addn2 = False, scaleP = 1.0)

    sim2.smartin.FWHM = res
    sim2.smartin.sample_res = res

    sim2.smartin.minwn = 1e4/lamax
    sim2.smartin.maxwn = 1e4/lamin 

    sim2.lblin.minwn = 1e4/lamax
    sim2.lblin.maxwn = 1e4/lamin

    if cia == "new":
        o2 = sim2.atmosphere.gases[2]
        o2.cia_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/o4_calc.cia'
    elif cia == "none":
        o2 = sim2.atmosphere.gases[2]
        o2.cia_file = None
    else:
        pass 

    sim2.gen_lblscripts()
    sim2.run_lblabc()
    sim2.write_smart(write_file = True)
    sim2.run_smart()

    sim2.open_outputs()
    wl2 = sim2.output.rad.lam
    flux2 = sim2.output.rad.pflux
    sflux2 = sim2.output.rad.sflux

    adj_flux2 = flux2/sflux2 * math.pi
    return(wl2, adj_flux2)

def integrate(lamin, lamax, atmos, cia):
    f = open("/gscratch/vsm/mwjl/projects/high_res/output/integrations_new.txt", "a")
    if atmos == 0:
        wl, flux = earth_like(lamin, lamax, 0.01, cia)
        wl_low, flux_low = earth_like(lamin, lamax,1, cia)
        tag = "el" #earth-like
    elif atmos == 1:
        wl, flux = ocean_loss(lamin, lamax, 0.01, cia)
        wl_low, flux_low = ocean_loss(lamin, lamax,1, cia)
        tag = "ol" # ocean loss
    else:
        wl, flux = ocean_outgassing(lamin, lamax, 0.01, cia)
        wl_low, flux_low = ocean_outgassing(lamin, lamax,1, cia)
        tag = "oo" # ocean outgassing
        
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
    while i < len(long_flux)- 25: 
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
    name = str(abs(adds)), str(lamin) + "to" + str(lamax), str(cia), str(tag)
    f = open("/gscratch/vsm/mwjl/projects/high_res/output/integrations_new.txt", "a")
    f.write(str(name) + "\n")

def output(lamin, lamax):
    integrate(lamin,lamax, 0, "new")
    integrate(lamin,lamax, 1, "new")
    integrate(lamin,lamax, 2, "new")

  #  integrate(lamin,lamax, 0, "none")
  #  integrate(lamin,lamax, 1, "none")
  #  integrate(lamin,lamax, 2, "none")

   # integrate(lamin,lamax, 0, "original")
   # integrate(lamin,lamax, 1, "original")
   # integrate(lamin,lamax, 2, "original")
 

if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="integ8",
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
         #output(0.78, 0.81)
         output(0.61, 0.65)
         output(0.67, 0.71)
         output(0.74, 0.78)
         output(1.24, 1.28)

    else:
        output(0.61, 0.65)
