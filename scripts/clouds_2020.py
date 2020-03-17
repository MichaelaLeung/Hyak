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
import scipy
from scipy.interpolate import interp1d
import scipy.integrate as integrate

def earth_like_hyak(lamin, lamax, res):
    name = 'prox'
    sim = smart.interface.Smart(tag = name)
    minwn = int(1e4/lamax)
    maxwn = int(1e4/lamin)
    smart_file = name + "_" + str(minwn) + "_" + str(maxwn) + "cm_toa.rad"
    try:
        f = open(smart_file)
        print("file exists")
        data = smart.readsmart.read_rad(smart_file)
        wl = data.lam
        flux = data.pflux
        sflux = data.sflux
        flux = flux/sflux
    except IOError:
        print("File does not exist")
        
        place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
        sim.set_run_in_place(place)
        
        sim.smartin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
        sim.lblin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
        sim.smartin.abs_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'

        infile = "/gscratch/vsm/mwjl/projects/high_res/inputs/profile_Earth_proxb_.pt_filtered"
        label = "Simulated Earth-like planet orbiting Proxima Centauri"
        sim.smartin.alb_file = "/gscratch/vsm/mwjl/projects/high_res/inputs/composite1_txt.txt"
        sim.set_planet_proxima_b()
        sim.load_atmosphere_from_pt(infile, addn2 = True)
        
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
        flux = (flux/sflux)
    return(wl, flux)

def earth(lamin, lamax, res):
    
    place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim = smart.interface.Smart(tag = "earth")
    sim.set_run_in_place(place)
    sim.smartin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.lblin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.smartin.abs_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'

    infile = "/gscratch/vsm/mwjl/projects/high_res/inputs/earth_avg.pt"
    label = "Simulated Earth-like planet orbiting Proxima Centauri"
    sim.smartin.alb_file = "/gscratch/vsm/mwjl/projects/high_res/inputs/composite1_txt.txt"
    sim.set_planet_proxima_b()
    sim.load_atmosphere_from_pt(infile, addn2 = True)
    
    o2 = sim.atmosphere.gases[3]
    o2.cia_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/o4_calc.cia'
    label = "Earth"

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

    sim.lblin.par_index = 7
    sim.smartin.out_level = 3
    sim.gen_lblscripts()
    sim.run_lblabc()
    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    sur_path_temp = str(sim.output.rad.path)
    sur_path = sur_path_temp[:-7] + "sur.rad"

    infile = sur_path

    # Convert each line to vector, compose array of vectors
    arrays = np.array([np.array(list(map(float, line.split()))) for line in open(infile)])

    # Flatten and reshape into rectangle grid
    arr = np.hstack(arrays).reshape((14, -1), order='F')

    # Parse columns
    lam   = arr[0,:]
    wno   = arr[1,:]
    solar = arr[2,:]
    dir_flux  = arr[3,:]
    diff_flux  = arr[4,:]

    total_flux = dir_flux + diff_flux 
    print(len(lam), len(total_flux))
    transmiss = total_flux / solar
    transmiss = transmiss / max(transmiss)
    return(lam, transmiss)


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


def cirrus(lamin, lamax, res):

    name = 'prox_cirrus'
    sim = smart.interface.Smart(tag = name)
    minwn = int(1e4/lamax)
    maxwn = int(1e4/lamin)
    smart_file = name + "_" + str(minwn) + "_" + str(maxwn) + "cm_toa.rad"
    try:
        f = open(smart_file)
        print("file exists")
        data = smart.readsmart.read_rad(smart_file)
        wl = data.lam
        flux = data.pflux
        sflux = data.sflux
        flux = flux/sflux
    except IOError:
        print("File does not exist")

        place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
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
                                              vert_xscale = 1.0e5,
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

        flux = flux/sflux

    return(wl, flux)

def strato(lamin, lamax, res):
    name = 'prox_strato'
    sim = smart.interface.Smart(tag = name)
    minwn = int(1e4/lamax)
    maxwn = int(1e4/lamin)
    smart_file = name + "_" + str(minwn) + "_" + str(maxwn) + "cm_toa.rad"
    try:
        f = open(smart_file)
        print("file exists")
        data = smart.readsmart.read_rad(smart_file)
        wl = data.lam
        flux = data.pflux
        sflux = data.sflux
        flux = flux/sflux
    except IOError:
        print("File does not exist")

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
                                              vert_xscale = 1.0e5,
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

        flux = flux/sflux

    return(wl, flux)

def ocean_outgassing(lamin, lamax, res):
    name = 'highw'
    sim = smart.interface.Smart(tag = name)
    minwn = int(1e4/lamax)
    maxwn = int(1e4/lamin)
    smart_file = name + "_" + str(minwn) + "_" + str(maxwn) + "cm_toa.rad"
    try:
        f = open(smart_file)
        print("file exists")
        data = smart.readsmart.read_rad(smart_file)
        wl = data.lam
        flux = data.pflux
        sflux = data.sflux
        flux = flux/sflux
    except IOError:
        print("File does not exist")

        place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'

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
        wl = sim2.output.rad.lam
        flux = sim2.output.rad.pflux
        sflux = sim2.output.rad.sflux

        flux = flux/sflux
    return(wl, flux)

def outgassing_cirrus(lamin, lamax, res):
    name = 'highw_cirrus'
    sim = smart.interface.Smart(tag = name)
    minwn = int(1e4/lamax)
    maxwn = int(1e4/lamin)
    smart_file = name + "_" + str(minwn) + "_" + str(maxwn) + "cm_toa.rad"
    try:
        f = open(smart_file)
        print("file exists")
        data = smart.readsmart.read_rad(smart_file)
        wl = data.lam
        flux = data.pflux
        sflux = data.sflux
        flux = flux/sflux
    except IOError:
        print("File does not exist")

        place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
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

        flux = flux/sflux

    return(wl, flux)

def outgassing_strato(lamin, lamax, res):
    name = 'highw_strato'
    sim = smart.interface.Smart(tag = name)
    minwn = int(1e4/lamax)
    maxwn = int(1e4/lamin)
    smart_file = name + "_" + str(minwn) + "_" + str(maxwn) + "cm_toa.rad"
    try:
        f = open(smart_file)
        print("file exists")
        data = smart.readsmart.read_rad(smart_file)
        wl = data.lam
        flux = data.pflux
        sflux = data.sflux
        flux = flux/sflux
    except IOError:
        print("File does not exist")
        place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
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

        flux = flux/sflux

    return(wl, flux)

def strato_noCIA(lamin, lamax, res):

    place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim = smart.interface.Smart(tag = "strato_noO4")
    sim.set_run_in_place(place)
    
    sim.smartin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.lblin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.smartin.abs_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'

    infile = "/gscratch/vsm/mwjl/projects/high_res/inputs/profile_Earth_proxb_.pt_filtered"
    label = "Simulated Earth-like planet orbiting Proxima Centauri"
    sim.smartin.alb_file = "/gscratch/vsm/mwjl/projects/high_res/inputs/composite1_txt.txt"
    sim.set_planet_proxima_b()
    sim.load_atmosphere_from_pt(infile, addn2 = True)
    
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

def cirrus_noCIA(lamin, lamax, res):

    place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim = smart.interface.Smart(tag = "cirrus_noO4")
    sim.set_run_in_place(place)
    
    sim.smartin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.lblin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.smartin.abs_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'

    infile = "/gscratch/vsm/mwjl/projects/high_res/inputs/profile_Earth_proxb_.pt_filtered"
    label = "Simulated Earth-like planet orbiting Proxima Centauri"
    sim.smartin.alb_file = "/gscratch/vsm/mwjl/projects/high_res/inputs/composite1_txt.txt"
    sim.set_planet_proxima_b()
    sim.load_atmosphere_from_pt(infile, addn2 = True)
    
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

def outgassing_strato_noCIA(lamin, lamax, res):

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

def outgassing_cirrus_noCIA(lamin, lamax, res):

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
    wl, flux = earth_like_hyak(lamin, lamax, 0.01)
    wl2, flux2 = cirrus(lamin, lamax, 0.01)
    wl3, flux3 = strato(lamin, lamax, 0.01)
    m, m_clouds = smart.utils.get_common_masks(wl, wl2)
    print("m,m_clouds", m, m_clouds)
    print("length of fluxes with m cloud mask applied", len(flux[m_clouds]), len(flux2[m_clouds]), len(flux3[m_clouds]))
    avg_flux = (0.5*flux[m_clouds]+0.25*flux2[m_clouds]+0.25*flux3[m_clouds])
    
    ocean_wl, ocean_flux = ocean_outgassing(lamin, lamax, 0.01)
    ocean_wl2, ocean_flux2 = outgassing_cirrus(lamin, lamax, 0.01)
    ocean_wl3, ocean_flux3 = outgassing_strato(lamin, lamax, 0.01)
    m, m_clouds = smart.utils.get_common_masks(ocean_wl, ocean_wl2)
    print("m, m_clouds (ocean)",m, m_clouds)
    print("length of ocean fluxes with m cloud mask applied",len(ocean_flux), len(ocean_flux2), len(ocean_flux3))
    avg_flux2 = (0.5*ocean_flux[m_clouds]+0.25*ocean_flux2[m_clouds]+0.25*ocean_flux3[m_clouds])

    fig, ax = plt.subplots(figsize = (10,10))
    ax.plot(wl[m_clouds], avg_flux[m_clouds], label = "1 bar Earth-Like")
    ax.plot(ocean_wl[m_clouds], avg_flux2[m_clouds], label = "10 bar Ocean Outgassing")
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

    ocean_wl, ocean_flux = ocean_outgassing(lamin, lamax, 0.01)
    ocean_wl2, ocean_flux2 = outgassing_cirrus(lamin, lamax, 0.01)
    ocean_wl3, ocean_flux3 = outgassing_strato(lamin, lamax, 0.01)
    m, m_clouds = smart.utils.get_common_masks(ocean_wl, ocean_wl2)
    avg_flux = (0.5*ocean_flux[m]+0.25*ocean_flux2[m]+0.25*ocean_flux3[m])
    
    noo4_wl, noo4_flux = ocean_outgassing(lamin, lamax, 0.01)
    noo4_wl2, noo4_flux2 = cirrus_noCIA(lamin, lamax, 0.01)
    noo4_wl3, noo4_flux3 = strato_noCIA(lamin, lamax,0.01)
    print(len(noo4_wl), len(noo4_wl2), len(noo4_wl3))
    m, m_clouds = smart.utils.get_common_masks(noo4_wl, noo4_wl3)
    m_clouds = m_clouds[:-1]
    avg_flux2 = (0.5*noo4_flux[m]+0.25*noo4_flux2[m]+0.25*noo4_flux3[m])

    
    fig, ax = plt.subplots(figsize = (10,10))
    print("oceanwl, avg flux", len(ocean_wl2), len(avg_flux))
    ax.plot(ocean_wl2[:len(avg_flux)], avg_flux[:len(ocean_wl2)], label = "10 bar Ocean Outgassing")
    ax.plot(noo4_wl[:len(avg_flux2)], avg_flux2[:len(avg_flux2)], label = "10 bar Ocean Loss, no O$_2$-O$_2$")
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
    
    ocean_wl, ocean_flux = ocean_outgassing(lamin, lamax,0.01)
    ocean_wl2, ocean_flux2 = outgassing_cirrus(lamin, lamax, 0.01)
    ocean_wl3, ocean_flux3 = outgassing_strato(lamin, lamax, 0.01)
    m, m_clouds = smart.utils.get_common_masks(ocean_wl, ocean_wl2)
    avg_flux = (0.5*ocean_flux[m]+0.25*ocean_flux2[m]+0.25*ocean_flux3[m])
    ax.plot(wl, avg_flux)
    ax.set_title("Simulated 10 bar oxygen ocean planet orbiting Proxima Centauri")
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$m)")

def cloud_weight_highw(lamin, lamax, res):
    ocean_wl, ocean_flux = ocean_outgassing(lamin, lamax, res)
    ocean_wl2, ocean_flux2 = outgassing_cirrus(lamin, lamax,res)
    ocean_wl3, ocean_flux3 = outgassing_strato(lamin, lamax, res)
    m, m_clouds = smart.utils.get_common_masks(ocean_wl, ocean_wl2)
    avg_flux = (0.5*ocean_flux[m_clouds]+0.25*ocean_flux2[m_clouds]+0.25*ocean_flux3[m_clouds])
    return(ocean_wl, avg_flux)

def cloud_weight_highw_noo4(lamin, lamax, res):
    noo4_wl, noo4_flux = ocean_outgassing(lamin, lamax, res)
    noo4_wl2, noo4_flux2 = cirrus_noCIA(lamin, lamax, res)
    noo4_wl3, noo4_flux3 = strato_noCIA(lamin, lamax, res)
    print(len(noo4_wl), len(noo4_wl2), len(noo4_wl3))
    m, m_clouds = smart.utils.get_common_masks(noo4_wl, noo4_wl2)
    avg_flux = (0.5*noo4_flux[m]+0.25*noo4_flux2[m]+0.25*noo4_flux3[m])
    return(noo4_wl, avg_flux)

def cloud_weight(lamin, lamax, res):
    ocean_wl, ocean_flux = earth_like_hyak(lamin, lamax, res)
    ocean_wl2, ocean_flux2 = cirrus(lamin, lamax,res)
    ocean_wl3, ocean_flux3 = strato(lamin, lamax, res)
    m, m_clouds = smart.utils.get_common_masks(ocean_wl, ocean_wl2)
    avg_flux = (0.5*ocean_flux[m]+0.25*ocean_flux2[m]+0.25*ocean_flux3[m])
    return(ocean_wl, avg_flux)

def cloud_weight_noo4(lamin, lamax, res):
    noo4_wl, noo4_flux = earth_like_hyak(lamin, lamax, res)
    noo4_wl2, noo4_flux2 = cirrus_noCIA(lamin, lamax, res)
    noo4_wl3, noo4_flux3 = strato_noCIA(lamin, lamax, res)
    print(len(noo4_wl), len(noo4_wl2), len(noo4_wl3))
    m, m_clouds = smart.utils.get_common_masks(noo4_wl, noo4_wl2)
    avg_flux = (0.5*noo4_flux[m]+0.25*noo4_flux2[m]+0.25*noo4_flux3[m])
    return(noo4_wl, avg_flux)
    
def high_pass(wl, flux, lamin, lamax, type):
    if type == 1:
        wl_low, flux_low = earth(lamin,lamax, 1)
    if type == 2: 
        wl_low, flux_low = cloud_weight(lamin,lamax, 1)

    long_flux = []
    for i in flux_low:
        j = 0
        while j < 101: 
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
    return(wl[:len(out)], out[:len(wl)])

def phase_calc(lamin,lamax):
    earth_wl, earth_flux = earth(lamin,lamax, 0.01)
    
    wl, flux = cloud_weight(lamin,lamax, 0.01)

    n_phase = 100
    phases = np.linspace(0,2*np.pi,n_phase)
    inclination = np.pi/2
    phi_90 = np.pi/2
    sma = 7500000
    c = 299792.458

    fluxes = np.outer(earth_flux, np.ones(n_phase))
    temp = np.arccos(-np.sin(inclination)*np.cos(phases))
    phase_function = ((np.sin(temp)+(np.pi-temp)*(np.cos(temp)))/(np.pi))
    #season = 0.55 #0-2PI march max
    season = -0.55 #july min?
    rv_bary = (29.8 * np.sin((11.2/365.)*phases + season))
    #rv_orb = (np.sqrt((G*m_prox)/(sma_m)) * np.sin(inclination)/1000 *np.sin(phases))
    rv_orb = (sma*2*np.pi)/967680* np.sin(inclination) *np.sin(phases)
    rv_sys = -21.7 * np.ones_like(phases)
    #rv_bary = (29.8 * np.sin((11.2/365.)*phases))
    rv = rv_sys + rv_orb - rv_bary
    inclination = np.pi/2.9
    obs_wl = np.outer(wl,(1+rv/c))

    i = 0
    out = []
    while i < n_phase:
        print("length of flux, obs_wl",len(flux), np.shape(obs_wl))
        high_pass_wl, high_pass_flux = high_pass(obs_wl[:len(flux),i],flux[:len(obs_wl)], lamin,lamax,2) 
        print("len high pass wl, flux", len(high_pass_wl), len(high_pass_flux))
        interp = scipy.interpolate.interp1d(high_pass_wl, high_pass_flux,fill_value = "extrapolate")
        out.append(interp(earth_wl) * earth_flux)
        i = i+1
    return(wl, out)

def integ_calc(lamin,lamax):
    wl, flux = cloud_weight(lamin, lamax, 0.01)
    wl, out = high_pass(wl, flux, lamin, lamax, 2)
    print("wl, out", np.shape(wl), np.shape(out))
    adds = integrate.trapz(out[:len(wl)], wl[:len(out)])

    print(adds)    
    name = str(abs(adds)) 
    f = open("/gscratch/vsm/mwjl/projects/high_res/scripts/integrations_" + str(lamin) + ".txt", "a")
    f.write(str(name) + "\n")
    f.close()
    
def master_plot():
    integ_calc(0.74, 0.78)
    integ_calc(0.67, 0.71)
    integ_calc(0.61, 0.65)
    integ_calc(1.25, 1.27)

def read_integ():
    print('starting read_integ')
    master_plot()
    files = "integrations_0.61.txt","integrations_0.67.txt", "integrations_0.74.txt", "integrations_1.25.txt"
    for name in files: 
       
        output = np.genfromtxt(name)
        phase = np.asarray(range(100))
        phase = phase* 2 * np.math.pi / 100
        print("output, phase", np.shape(output), np.shape(phase))
        phase.astype(np.float)
#        fig, ax = plt.subplots(figsize = (12,12))
#        ax.plot(phase, output)
#        ax.set_title("Integration Metric over Phase")
#        ax.set_ylabel("Integration Metric")
#        ax.set_xlabel("Phase")
#        out_name = name[:-4] + ".png"
#        fig.savefig(out_name, bbox_inches = 'tight')
        integ = integrate.trapz(output,phase) 
        f = open("/gscratch/vsm/mwjl/projects/high_res/scripts/integrations_fin_w.txt", "a")
        f.write(str(integ) + '\n')
        f.close()
    


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
       # plotting(0.61,0.65, "Gamma band (0.63) Ocean Outgassing")
      #  plotting(0.67,0.71, "Oxygen B band (0.69) Ocean Outgassing")
      #  plotting(0.74,0.78,"Oxygen A band (0.76) Ocean Outgassing")
      #  plotting(1.25,1.29,"1.27 Ocean Outgassing")
        
#        integ_calc(0.61, 0.65)
#        integ_calc(0.67, 0.71)
#        integ_calc(0.74, 0.78)
#        integ_calc(1.25,1.29)

#        plotting_noO4(0.61,0.65, "Gamma band (0.63) Ocean Outgassing")
#        plotting_noO4(0.67,0.71, "Oxygen B band (0.69) Ocean Outgassing")
#        plotting_noO4(0.74,0.78,"Oxygen A band (0.76) Ocean Outgassing")
#        plotting_noO4(1.25,1.29,"1.27 Ocean Outgassing")

        long()
#        master_plot()
#        read_integ()
        
    else:
        plotting(0.61,0.645,1,"Gamma band (0.63) Ocean Outgassing")

