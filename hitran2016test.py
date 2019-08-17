#!/usr/bin/python
import numpy as np
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.collections import LineCollection
from astropy.io import fits 
import smart
import sys, os
import datetime
matplotlib.rcParams['text.usetex'] = False
import random
import platform
lamin = 0.7
lamax = 0.8
res = 1/(10*lamin)

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
o2.cia_file = '/gscratch/vsm/mwjl/projects/high_res/smart_output/cia_adj_calc.cia'

label = "Earth-Like"
sim.set_planet_proxima_b()
sim.set_star_proxima()

sim.set_run_in_place() 
sim.set_executables_automatically()

sim.lblin.par_file = '/gscratch/vsm/alinc/fixed_input/HITRAN2016.par' #/gscratch/vsm/alinc/fixed_input/
sim.lblin.hitran_tag = 'hitran2016'
sim.lblin.fundamntl_file = '/gscratch/vsm/alinc/fixed_input/fundamntl2016.dat'
sim.lblin.lblabc_exe = '/gscratch/vsm/alinc/exec/lblabc_2016'

sim.smartin.sza = 57

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

sim.open_outputs()
wl = sim.output.rad.lam
flux = sim.output.rad.pflux
sflux = sim.output.rad.sflux
adj_flux = flux/sflux
