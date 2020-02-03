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
from scipy.interpolate import interp1d

def earth_like_hyak(lamin, lamax):
    
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
    o2.cia_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/o4_calc.cia'
    label = "Earth-Like"
    sim.set_planet_proxima_b()
    sim.set_star_proxima()

    sim.set_executables_automatically()

    sim.lblin.par_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/HITRAN2019.par' #/gscratch/vsm/alinc/fixed_input/
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
    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux
    adj_flux = (flux/sflux)
    return(wl, adj_flux)

def cloud_frac():
    infile2 = "/gscratch/vsm/mwjl/projects/high_res/inputs/10bar_O2_wet.pt_filtered.pt"
    data = np.genfromtxt(infile2)
    data = data[56:61, 0:2]
    press = data[:,0]
    press = np.log(press)
    temp = data[:,1]
    interpd = interp1d(temp, press)
    strato_temp = 281.1
    cirrus_temp = 227.5
    strato_out = interpd(strato_temp)
    cirrus_out = interpd(cirrus_temp)
    strato = np.exp(strato_out)
    cirrus = np.exp(cirrus_out)
    strato = strato / 10**5
    cirrus = cirrus / 10**5

    strato_bar = 0.847
    cirrus_bar = 0.331

    f_cirrus = cirrus / cirrus_bar * 10**5
    f_strato = strato / strato_bar * 10**5

    return(f_cirrus, f_strato)


def cirrus(lamin, lamax):
    res = 1/(10*lamin)

    place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim = smart.interface.Smart(tag = "prox_cirrus")
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

    sim.lblin.par_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/HITRAN2019.par'
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

    f_cirrus, f_strato = cloud_frac()
    
    # Create a cirrus cloud mie scattering aerosol mode
    mie_cirrus = smart.interface.MieMode(mie_file = os.path.join(smart.interface.CLDMIEDIR, "baum_cirrus_de100.mie"),
                                         mie_skip = 1,
                                         mie_lines =
                                         '1,4,5,3',
                                         iang_smart = 2)

    # Create an optical depth profile
    tau_cirrus = smart.interface.CloudTau(vert_file = os.path.join(smart.interface.CLDMIEDIR, "cld_tau.dat"),
                                          vert_ref_wno = 1e4/lamax,
                                          vert_skip = 4,
                                          vert_coord = 1,
                                          vert_xscale = f_cirrus,
                                          vert_yscale = 2.0)

    # Create an Aerosol object with our cirrus mie scattering and optical depths
    cirrus = smart.interface.Aerosols(miemodes=[mie_cirrus],
                                      mietau=[tau_cirrus])

    sim.aerosols = cirrus

    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux

    adj_flux = flux/sflux

    return(wl, adj_flux)

def strato(lamin, lamax):
    res = 1/(10*lamin)

    place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim = smart.interface.Smart(tag = "prox_strato")
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

    sim.lblin.par_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/HITRAN2019.par'
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

    f_cirrus, f_strato = cloud_frac()

    # Create a stratocumulus cloud mie scattering aerosol mode
    mie_strato = smart.interface.MieMode(mie_file = os.path.join(smart.interface.CLDMIEDIR, "strato_cum.mie"),
                                         mie_skip = 19,
                                         mie_lines = '1,7,8,11',
                                         iang_smart = 1,
                                         mom_skip = 17)

    # Create an optical depth profile
    tau_strato = smart.interface.CloudTau(vert_file = os.path.join(smart.interface.CLDMIEDIR, "cld_tau.dat"),
                                          vert_ref_wno = 1e4/lamax,
                                          vert_skip = 28,
                                          vert_coord = 1,
                                          vert_xscale = f_strato,
                                          vert_yscale = 1.0)

    # Create an Aerosol object with our stratocumulus mie scattering and optical depths
    strato = smart.interface.Aerosols(miemodes=[mie_strato],
                                      mietau=[tau_strato])

    sim.aerosols = strato


    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux

    adj_flux = flux/sflux

    return(wl, adj_flux)

def ocean_outgassing(lamin, lamax):
    res = 1/(10*lamin)
    place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'

    sim2 = smart.interface.Smart(tag = "highw")
    sim2.set_run_in_place(place)
    sim2.smartin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim2.lblin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim2.smartin.abs_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    infile2 = "/gscratch/vsm/mwjl/projects/high_res/inputs/10bar_O2_wet.pt_filtered.pt"
    label = "Ocean Outgassing"
    sim2.smartin.alb_file = "/gscratch/vsm/mwjl/projects/high_res/inputs/earth_noveg_highw.alb"
    sim2.set_planet_proxima_b()
    sim2.set_star_proxima()

    sim2.set_run_in_place() 
    sim2.set_executables_automatically()

    sim2.lblin.par_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/HITRAN2019.par' #/gscratch/vsm/alinc/fixed_input/
    sim2.lblin.hitran_tag = 'hitran2016'
    sim2.lblin.fundamntl_file = '/gscratch/vsm/alinc/fixed_input/fundamntl2016.dat'
    sim2.lblin.lblabc_exe = '/gscratch/vsm/alinc/exec/lblabc_2016'
    sim2.lblin.par_index = 7

    sim2.smartin.iraylei = [4]
    sim2.smartin.vraylei = [1]

    sim2.smartin.sza = 57
    sim2.load_atmosphere_from_pt(infile2, addn2 = False, scaleP = 1.0)

    sim2.smartin.FWHM = res
    sim2.smartin.sample_res = res

    sim2.smartin.minwn = 1e4/lamax
    sim2.smartin.maxwn = 1e4/lamin 

    sim2.lblin.minwn = 1e4/lamax
    sim2.lblin.maxwn = 1e4/lamin 


    o2 = sim2.atmosphere.gases[2]
    o2.cia_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/o4_calc.cia'

    sim2.gen_lblscripts()
    sim2.run_lblabc()
    sim2.write_smart(write_file = True)
    sim2.run_smart()

    sim2.open_outputs()
    wl2 = sim2.output.rad.lam
    flux2 = sim2.output.rad.pflux
    sflux2 = sim2.output.rad.sflux

    adj_flux2 = flux2/sflux2
    return(wl2, adj_flux2)

def outgassing_cirrus(lamin, lamax):
    res = 1/(10*lamin)

    place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim = smart.interface.Smart(tag = "highw_cirrus")
    sim.set_run_in_place(place)
    
    sim.smartin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.lblin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.smartin.abs_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'

    infile = "/gscratch/vsm/mwjl/projects/high_res/inputs/10bar_O2_wet.pt_filtered.pt"
    label = "Ocean Outgassing"
    sim.smartin.alb_file = "/gscratch/vsm/mwjl/projects/high_res/inputs/earth_noveg_highw.alb"
    sim.set_planet_proxima_b()
    sim.load_atmosphere_from_pt(infile, addn2 = False)
    
    o2 = sim.atmosphere.gases[3]
    o2.cia_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/o4_calc.cia'
    label = "Earth-Like"
    sim.set_planet_proxima_b()
    sim.set_star_proxima()

    sim.set_executables_automatically()

    sim.lblin.par_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/HITRAN2019.par'
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

    f_cirrus, f_strato = cloud_frac()
    
    # Create a cirrus cloud mie scattering aerosol mode
    mie_cirrus = smart.interface.MieMode(mie_file = os.path.join(smart.interface.CLDMIEDIR, "baum_cirrus_de100.mie"),
                                         mie_skip = 1,
                                         mie_lines =
                                         '1,4,5,3',
                                         iang_smart = 2)

    # Create an optical depth profile
    tau_cirrus = smart.interface.CloudTau(vert_file = os.path.join(smart.interface.CLDMIEDIR, "cld_tau.dat"),
                                          vert_ref_wno = 1e4/lamax,
                                          vert_skip = 4,
                                          vert_coord = 1,
                                          vert_xscale = f_cirrus,
                                          vert_yscale = 2.0)

    # Create an Aerosol object with our cirrus mie scattering and optical depths
    cirrus = smart.interface.Aerosols(miemodes=[mie_cirrus],
                                      mietau=[tau_cirrus])

    sim.aerosols = cirrus

    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux

    adj_flux = flux/sflux

    return(wl, adj_flux)

def outgassing_strato(lamin, lamax):
    res = 1/(10*lamin)

    place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim = smart.interface.Smart(tag = "highw_strato")
    sim.set_run_in_place(place)
    
    sim.smartin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.lblin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.smartin.abs_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'

    infile = "/gscratch/vsm/mwjl/projects/high_res/inputs/10bar_O2_wet.pt_filtered.pt"
    label = "Ocean Outgassing"
    sim.smartin.alb_file = "/gscratch/vsm/mwjl/projects/high_res/inputs/earth_noveg_highw.alb"
    sim.set_planet_proxima_b()
    sim.load_atmosphere_from_pt(infile, addn2 = False)
    
    o2 = sim.atmosphere.gases[3]
    o2.cia_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/o4_calc.cia'
    label = "Earth-Like"
    sim.set_planet_proxima_b()
    sim.set_star_proxima()

    sim.set_executables_automatically()

    sim.lblin.par_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/HITRAN2019.par'
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

    f_cirrus, f_strato = cloud_frac()
    
     # Create a stratocumulus cloud mie scattering aerosol mode
    mie_strato = smart.interface.MieMode(mie_file = os.path.join(smart.interface.CLDMIEDIR, "strato_cum.mie"),
                                         mie_skip = 19,
                                         mie_lines = '1,7,8,11',
                                         iang_smart = 1,
                                         mom_skip = 17)

    # Create an optical depth profile
    tau_strato = smart.interface.CloudTau(vert_file = os.path.join(smart.interface.CLDMIEDIR, "cld_tau.dat"),
                                          vert_ref_wno = 1e4/lamax,
                                          vert_skip = 28,
                                          vert_coord = 1,
                                          vert_xscale = f_strato,
                                          vert_yscale = 1.0)

    # Create an Aerosol object with our stratocumulus mie scattering and optical depths
    strato = smart.interface.Aerosols(miemodes=[mie_strato],
                                      mietau=[tau_strato])

    sim.aerosols = strato

    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux

    adj_flux = flux/sflux

    return(wl, adj_flux)

def strato_noCIA(lamin, lamax):
    res = 1/(10*lamin)

    place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim = smart.interface.Smart(tag = "highw_strato_noO4")
    sim.set_run_in_place(place)
    
    sim.smartin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.lblin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.smartin.abs_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'

    infile = "/gscratch/vsm/mwjl/projects/high_res/inputs/10bar_O2_wet.pt_filtered.pt"
    label = "Ocean Outgassing"
    sim.smartin.alb_file = "/gscratch/vsm/mwjl/projects/high_res/inputs/earth_noveg_highw.alb"
    sim.set_planet_proxima_b()
    sim.load_atmosphere_from_pt(infile, addn2 = False)
    
    o2 = sim.atmosphere.gases[3]
    o2.cia_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/o4_calc.cia'
    label = "Earth-Like"
    sim.set_planet_proxima_b()
    sim.set_star_proxima()

    sim.set_executables_automatically()

    sim.lblin.par_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/HITRAN2019.par'
    sim.lblin.hitran_tag = 'hitran2016'
    sim.lblin.fundamntl_file = '/gscratch/vsm/alinc/fixed_input/fundamntl2016.dat'
    sim.lblin.lblabc_exe = '/gscratch/vsm/alinc/exec/lblabc_2016'
    sim.lblin.par_index = 7


    sim.smartin.sza = 57

    o2 = sim.atmosphere.gases[2]
    o2.cia_file = None


    sim.smartin.FWHM = res
    sim.smartin.sample_res = res

    sim.smartin.minwn = 1e4/lamax
    sim.smartin.maxwn = 1e4/lamin 

    sim.lblin.minwn = 1e4/lamax
    sim.lblin.maxwn = 1e4/lamin 


    sim.gen_lblscripts()
    sim.run_lblabc()

    f_cirrus, f_strato = cloud_frac()

    
    # Create a cirrus cloud mie scattering aerosol mode
    mie_cirrus = smart.interface.MieMode(mie_file = os.path.join(smart.interface.CLDMIEDIR, "baum_cirrus_de100.mie"),
                                         mie_skip = 1,
                                         mie_lines =
                                         '1,4,5,3',
                                         iang_smart = 2)

    # Create an optical depth profile
    tau_cirrus = smart.interface.CloudTau(vert_file = os.path.join(smart.interface.CLDMIEDIR, "cld_tau.dat"),
                                          vert_ref_wno = 1e4/lamax,
                                          vert_skip = 4,
                                          vert_coord = 1,
                                          vert_xscale = f_cirrus,
                                          vert_yscale = 2.0)

    # Create an Aerosol object with our cirrus mie scattering and optical depths
    cirrus = smart.interface.Aerosols(miemodes=[mie_cirrus],
                                      mietau=[tau_cirrus])

    sim.aerosols = cirrus

    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux

    adj_flux = flux/sflux

    return(wl, adj_flux)

def cirrus_noCIA(lamin, lamax):
    res = 1/(10*lamin)

    place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim = smart.interface.Smart(tag = "highw_cirrus_noO4")
    sim.set_run_in_place(place)
    
    sim.smartin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.lblin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.smartin.abs_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'

    infile = "/gscratch/vsm/mwjl/projects/high_res/inputs/10bar_O2_wet.pt_filtered.pt"
    label = "Ocean Outgassing"
    sim.smartin.alb_file = "/gscratch/vsm/mwjl/projects/high_res/inputs/earth_noveg_highw.alb"
    sim.set_planet_proxima_b()
    sim.load_atmosphere_from_pt(infile, addn2 = False)
    
    o2 = sim.atmosphere.gases[3]
    o2.cia_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/o4_calc.cia'
    label = "Earth-Like"
    sim.set_planet_proxima_b()
    sim.set_star_proxima()

    sim.set_executables_automatically()

    sim.lblin.par_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/HITRAN2019.par'
    sim.lblin.hitran_tag = 'hitran2016'
    sim.lblin.fundamntl_file = '/gscratch/vsm/alinc/fixed_input/fundamntl2016.dat'
    sim.lblin.lblabc_exe = '/gscratch/vsm/alinc/exec/lblabc_2016'
    sim.lblin.par_index = 7


    sim.smartin.sza = 57

    o2 = sim.atmosphere.gases[2]
    o2.cia_file = None


    sim.smartin.FWHM = res
    sim.smartin.sample_res = res

    sim.smartin.minwn = 1e4/lamax
    sim.smartin.maxwn = 1e4/lamin 

    sim.lblin.minwn = 1e4/lamax
    sim.lblin.maxwn = 1e4/lamin 


    sim.gen_lblscripts()
    sim.run_lblabc()

    f_cirrus, f_strato = cloud_frac()

    
    # Create a cirrus cloud mie scattering aerosol mode
    mie_cirrus = smart.interface.MieMode(mie_file = os.path.join(smart.interface.CLDMIEDIR, "baum_cirrus_de100.mie"),
                                         mie_skip = 1,
                                         mie_lines =
                                         '1,4,5,3',
                                         iang_smart = 2)

    # Create an optical depth profile
    tau_cirrus = smart.interface.CloudTau(vert_file = os.path.join(smart.interface.CLDMIEDIR, "cld_tau.dat"),
                                          vert_ref_wno = 1e4/lamax,
                                          vert_skip = 4,
                                          vert_coord = 1,
                                          vert_xscale = f_cirrus,
                                          vert_yscale = 2.0)

    # Create an Aerosol object with our cirrus mie scattering and optical depths
    cirrus = smart.interface.Aerosols(miemodes=[mie_cirrus],
                                      mietau=[tau_cirrus])

    sim.aerosols = cirrus

    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux

    adj_flux = flux/sflux

    return(wl, adj_flux)



def plotting(lamin, lamax, title):
    matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})
    matplotlib.rcParams['font.size'] = 25.0
    matplotlib.rc('text', usetex=False)
    plt.switch_backend('agg')
    fig_name = int(100*(float(lamin) + float(lamax))/2)
    wl, flux = earth_like_hyak(lamin, lamax)
    wl2, flux2 = cirrus(lamin, lamax)
    wl3, flux3 = strato(lamin, lamax)
    print(len(flux), len(flux2), len(flux3))
    avg_flux = (0.5*flux+0.25*flux2+0.25*flux3)
    
    ocean_wl, ocean_flux = ocean_outgassing(lamin, lamax)
    ocean_wl2, ocean_flux2 = outgassing_cirrus(lamin, lamax)
    ocean_wl3, ocean_flux3 = outgassing_strato(lamin, lamax)
    print(len(ocean_flux), len(ocean_flux2), len(ocean_flux3))
    avg_flux2 = (0.5*ocean_flux+0.25*ocean_flux2+0.25*ocean_flux3)

    fig, ax = plt.subplots(figsize = (10,10))
    ax.plot(wl, avg_flux, label = "1 bar Earth-Like")
    ax.plot(ocean_wl, avg_flux2, label = "10 bar Ocean Outgassing")
    ax.set_title(title)
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$m)")
    if lamin == 0.61:
        ax.legend()
    fig.savefig("/gscratch/vsm/mwjl/projects/high_res/plots/" + str(fig_name) +  "new_CIA_ocean_clouds.png", bbox_inches = "tight")


def plotting_noO4(lamin, lamax,  title):
    matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})
    matplotlib.rcParams['font.size'] = 25.0
    matplotlib.rc('text', usetex=False)
    plt.switch_backend('agg')
    fig_name = int(100*(float(lamin) + float(lamax))/2)
    ocean_wl, ocean_flux = ocean_outgassing(lamin, lamax)
    ocean_wl2, ocean_flux2 = outgassing_cirrus(lamin, lamax)
    ocean_wl3, ocean_flux3 = outgassing_strato(lamin, lamax)
    avg_flux = (0.5*ocean_flux+0.25*ocean_flux2+0.25*ocean_flux3)
    
    noo4_wl, noo4_flux = ocean_outgassing(lamin, lamax)
    noo4_wl2, noo4_flux2 = cirrus_noCIA(lamin, lamax)
    noo4_wl3, noo4_flux3 = strato_noCIA(lamin, lamax)
    print(len(noo4_wl), len(noo4_wl2), len(noo4_wl3))
    avg_flux2 = (0.5*noo4_flux+0.25*noo4_flux2+0.25*noo4_flux3)

    
    fig, ax = plt.subplots(figsize = (10,10))
    ax.plot(ocean_wl, avg_flux, label = "10 bar Ocean Outgassing")
    ax.plot(noo4_wl, avg_flux2, label = "10 bar Ocean Loss, no O$_2$-O$_2$")
    ax.set_title(title)
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$m)")
    if lamin == 0.61:
        ax.legend()
    fig.savefig("/gscratch/vsm/mwjl/projects/high_res/plots/" + str(fig_name) +  "new_CIA_ocean_clouds_noO4.png", bbox_inches = "tight")


def long():
    lamin = 0.4
    lamax = 2.0
    
    fig, ax = plt.subplots(figsize = (30,10))
    
    ocean_wl, ocean_flux = ocean_outgassing(lamin, lamax)
    ocean_wl2, ocean_flux2 = outgassing_cirrus(lamin, lamax)
    ocean_wl3, ocean_flux3 = outgassing_strato(lamin, lamax)
    avg_flux = (0.5*ocean_flux+0.25*ocean_flux2+0.25*ocean_flux3)

    ax.plot(wl, avg_flux)
    ax.set_title("Simulated 10 bar oxygen ocean planet orbiting Proxima Centauri")
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$m)")
    


if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="nor_cld",
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
        plotting(0.61,0.65, "Gamma band (0.63) Ocean Outgassing")
        plotting(0.67,0.71, "Oxygen B band (0.69) Ocean Outgassing")
        plotting(0.74,0.78,"Oxygen A band (0.76) Ocean Outgassing")
        plotting(1.25,1.29,"1.27 Ocean Outgassing")

        plotting_noO4(0.61,0.65, "Gamma band (0.63) Ocean Outgassing")
        plotting_noO4(0.67,0.71, "Oxygen B band (0.69) Ocean Outgassing")
        plotting_nO4(0.74,0.78,"Oxygen A band (0.76) Ocean Outgassing")
        plotting_noO4(1.25,1.29,"1.27 Ocean Outgassing")
        long()
    else:
        plotting(0.61,0.645,1,"Gamma band (0.63) Ocean Outgassing")

