import numpy as np
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from astropy.io import fits
import smart
import sys, os
import datetime
matplotlib.rcParams['text.usetex'] = False
matplotlib.rcParams['legend.loc'] = 'lower center'
import random
import math
import platform 
import scipy
print('import finished') 
import scipy.integrate as integrate

import matplotlib; matplotlib.use('agg')
matplotlib.rcParams['text.usetex'] = False
import scipy
from scipy.interpolate import interp1d
import scipy.integrate as integrate

def earth_like_hyak(lamin, lamax, res):
    name = '/gscratch/vsm/mwjl/projects/high_res/smart_output/'+'prox' + str(res)
    sim = smart.interface.Smart(tag = 'prox'+str(res))
    minwn = int(1e4/lamax)
    maxwn = int(1e4/lamin)
    smart_file = name + "_" + str(minwn) + "_" + str(maxwn) + "cm_toa.rad"
    print(smart_file)
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
    else:
        data = smart.readsmart.read_rad(smart_file)
        wl = data.lam
        flux = data.pflux
        sflux = data.sflux
        return(wl, flux, sflux)
        f = open('/gscratch/vsm/mwjl/projects/high_res/smart_output/'+smart_file)
        print('/gscratch/vsm/mwjl/projects/high_res/smart_output/'+smart_file)
        return(wl, flux, sflux) 

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

    name = '/gscratch/vsm/mwjl/projects/high_res/smart_output/'+'prox_cirrus'+str(res)
    sim = smart.interface.Smart(tag = 'prox_cirrus'+str(res))
    minwn = int(1e4/lamax)
    maxwn = int(1e4/lamin)
    smart_file = name + "_" + str(minwn) + "_" + str(maxwn) + "cm_toa.rad"
    print(smart_file)
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
        return(wl, flux, sflux)

    else:
        data = smart.readsmart.read_rad(smart_file)
        wl = data.lam
        flux = data.pflux
        sflux = data.sflux
        flux = flux/sflux
        return(wl, flux, sflux)
        f = open('/gscratch/vsm/mwjl/projects/high_res/smart_output/'+smart_file)
        print('/gscratch/vsm/mwjl/projects/high_res/smart_output/'+smart_file)
        return(wl, flux, sflux)


    return(wl, flux)

def strato(lamin, lamax, res):
    name = '/gscratch/vsm/mwjl/projects/high_res/smart_output/'+'prox_strato'+str(res)
    sim = smart.interface.Smart(tag = 'prox_strato'+str(res))
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
        return(wl, flux, sflux)
    else:
        data = smart.readsmart.read_rad(smart_file)
        wl = data.lam
        flux = data.pflux
        sflux = data.sflux
        return(wl, flux, sflux)
        f = open('/gscratch/vsm/mwjl/projects/high_res/smart_output/'+smart_file)
        print('/gscratch/vsm/mwjl/projects/high_res/smart_output/'+smart_file)
        return(wl, flux, sflux)


    return(wl, flux)

def ocean_outgassing(lamin, lamax, res):
    name = '/gscratch/vsm/mwjl/projects/high_res/smart_output/'+'highw'+str(res)
    sim2 = smart.interface.Smart(tag = 'highw'+str(res))
    minwn = int(1e4/lamax)
    maxwn = int(1e4/lamin)
    smart_file = name + "_" + str(minwn) + "_" + str(maxwn) + "cm_toa.rad"
    try:
        f = open("/gscratch/vsm/mwjl/projects/high_res/smart_output/"+smart_file)
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

    return(wl, flux, sflux)

def outgassing_cirrus(lamin, lamax, res):
    name = '/gscratch/vsm/mwjl/projects/high_res/smart_output/'+'highw_cirrus'+str(res)
    sim = smart.interface.Smart(tag = 'highw_cirrus'+str(res))
    minwn = int(1e4/lamax)
    maxwn = int(1e4/lamin)
    smart_file = name + "_" + str(minwn) + "_" + str(maxwn) + "cm_toa.rad"
    try:
        f = open("/gscratch/vsm/mwjl/projects/high_res/smart_output/" + smart_file)
        print("file exists")
        data = smart.readsmart.read_rad(smart_file)
        wl = data.lam
        flux = data.pflux
        sflux = data.sflux
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
        sim.load_atmosphere_from_pt(infile, addn2 = False, scaleP = 1.0)
        
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

    return(wl, flux, sflux)

def outgassing_strato(lamin, lamax, res):
    name = '/gscratch/vsm/mwjl/projects/high_res/smart_output/'+'highw_strato'+str(res)
    sim = smart.interface.Smart(tag = 'highw_strato'+str(res))
    minwn = int(1e4/lamax)
    maxwn = int(1e4/lamin)
    smart_file = name + "_" + str(minwn) + "_" + str(maxwn) + "cm_toa.rad"
    try:
        f = open("/gscratch/vsm/mwjl/projects/high_res/smart_output/" + smart_file)
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
        sim.load_atmosphere_from_pt(infile, addn2 = False, scaleP = 1.0)
        
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

    return(wl, flux, sflux)

def cloud_weight_highw(lamin, lamax, res):
    ocean_wl, ocean_flux, ocean_sflux = ocean_outgassing(lamin, lamax, res)
    ocean_wl2, ocean_flux2, ocean_sflux2 = outgassing_cirrus(lamin, lamax,res)
    ocean_wl3, ocean_flux3, ocean_sflux3 = outgassing_strato(lamin, lamax, res)
    m, m_clouds = smart.utils.get_common_masks(ocean_wl, ocean_wl3)
    print(len(ocean_wl), len(ocean_wl2), len(ocean_wl3))
    avg_flux = (0.5*ocean_flux[m_clouds]+0.25*ocean_flux2[m_clouds]+0.25*ocean_flux3[m_clouds])
    avg_sflux = (0.5*ocean_sflux[m_clouds]+0.25*ocean_sflux2[m_clouds]+0.25*ocean_sflux3[m_clouds])
    return(ocean_wl, avg_flux, avg_sflux)

def cloud_weight(lamin, lamax, res):
    ocean_wl, ocean_flux, ocean_sflux= earth_like_hyak(lamin, lamax, res)
    ocean_wl2, ocean_flux2, ocean_sflux = cirrus(lamin, lamax,res)
    ocean_wl3, ocean_flux3, ocean_sflux = strato(lamin, lamax, res)
    m, m_clouds = smart.utils.get_common_masks(ocean_wl, ocean_wl2)
    avg_flux = (0.5*ocean_flux[m]+0.25*ocean_flux2[m]+0.25*ocean_flux3[m])
    avg_sflux = (0.5*ocean_flux[m]+0.25*ocean_flux2[m]+0.25*ocean_flux3[m])
    return(ocean_wl, avg_flux, avg_sflux)

def plotting(atmos):
    lamin = 0.5
    lamax = 2.0
    res = 0.01
    sim = smart.interface.Smart(tag = 'highw_strato')
    sim.set_planet_proxima_b()    
    r_km = 149598000 * sim.smartin.r_AU

    
    fig, ax = plt.subplots(figsize = (30, 10))
    if atmos == 'prox': 
        wl, flux, sflux = cloud_weight(lamin, lamax, res)
        label = "Simulated Earth-like planet orbiting Proxima Centauri"
    elif atmos == 'highw':
        wl, flux, sflux = cloud_weight_highw(lamin, lamax, res)
        label = "Simulated cloudy ocean outgassing planet orbiting Proxima Centauri"

    temp = flux / sflux 
    fpfs = flux/sflux * (sim.smartin.radius/r_km)**2

    ax.plot(wl[:len(flux)], temp[:len(wl)])
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$m)")
    ax2 = ax.twinx()
    ax2.set_ylabel("Planet-to-star contrast ratio")     
    ax2.plot(wl[:len(flux)], fpfs[:len(wl)])
    ax.set_title(label)
    ax.set_xlim(0.5,2)
    ax.axvspan(0.61, 0.65, alpha=0.5, color='0.85')
    ax.axvspan(0.67, 0.71, alpha=0.5, color='0.85')
    ax.axvspan(0.74, 0.78, alpha=0.5, color='0.85')
    ax.axvspan(1.25, 1.29, alpha=0.5, color='0.85')
    fig.savefig("/gscratch/vsm/mwjl/projects/high_res/plots/" + str(atmos) + "_newCIA.png", bbox_inches = 'tight')

if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="co_test",
                               subname="submit.csh",
                               workdir = "",
                               nodes = 1,
                               mem = "500G",
                               walltime = "24:00:00",
                               ntasks = 28,
                               account = "vsm",
                               submit = True,
                               rm_after_submit = True)
    elif platform.node().startswith("n"):
        # On a mox compute node: ready to run
        plotting("prox")
#        plotting("highd") 
        plotting("highw")
#    else:
        # Presumably, on a regular computer: ready to run
        plotting("prox")      









