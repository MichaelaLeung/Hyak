#!/usr/bin/python

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
import random
import platform

lamin = 0.7
lamax = 0.8
res = 1/(10*lamin)

sim = smart.interface.Smart(tag = "prox")
infile = "profile_Earth_proxb_.pt_filtered"
label = "Simulated Earth-like planet orbiting Proxima Centauri"
sim.smartin.alb_file = "composite1_txt.txt"
sim.set_planet_proxima_b()
sim.load_atmosphere_from_pt(infile, addn2 = False)
    
o2 = sim.atmosphere.gases[3]
o2.cia_file = 'cia_adj_calc.cia'
infile = "profile_Earth_proxb_.pt_filtered"
label = "Earth-Like"
sim.smartin.alb_file = "composite1_txt.txt"
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
    
mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})
mpl.rcParams['font.size'] = 25.0
mpl.rc('text', usetex=False)
plt.switch_backend('agg')
fig_name = int(100*(float(lamin) + float(lamax))/2)
fig, ax = plt.subplots(figsize = (10,10))
ax.plot(wl, adj_flux, label = "1 bar Earth-Like")
ax.set_ylabel("Reflectance")
ax.set_xlabel("Wavelength ($\mu$ m)")
ax.legend()
fig.savefig(str(fig_name) +  "2016_test.png", bbox_inches = "tight")
