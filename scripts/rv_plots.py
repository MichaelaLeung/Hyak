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
import imageio
import platform 
import scipy
from scipy import interpolate
print('import finished') 
import scipy.integrate as integrate

from clouds_2020 import cloud_weight_highw

high_res = 0.1

low_res = 10

n_phase = 10

def run_prox(lamin, lamax, res):

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
    sim.load_atmosphere_from_pt(infile, addn2 = True)
    
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
    sim.smartin.out_level = 1
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
    adj_flux = math.pi * (flux/sflux)
    return(wl, flux)

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

    adj_flux = flux/sflux

    return(wl, flux)

def run_earth(lamin, lamax, res):
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

def clouds_out(lamin, lamax, res):
    wl, flux = run_prox(lamin, lamax, res)
    wl2, flux2 = clouds(lamin, lamax, 0, res)
    wl3, flux3 = clouds(lamin, lamax, 1, res)
    avg_flux = (0.5*flux[:min(len(flux), len(flux2), len(flux3))]+0.25*flux2[:min(len(flux), len(flux2), len(flux3))]+0.25*flux3[:min(len(flux), len(flux2), len(flux3))])
    return(wl, avg_flux)

def ocean_loss(lamin, lamax, res):
    place = '/gscratch/vsm/mwjl/projects/high_res/smart_output'

    sim = smart.interface.Smart(tag = "highd")
    sim.set_run_in_place(place)
    sim.smartin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.lblin.out_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    sim.smartin.abs_dir = '/gscratch/vsm/mwjl/projects/high_res/smart_output'
    
    infile = "/gscratch/vsm/mwjl/projects/high_res/inputs/10bar_O2_dry.pt_filtered.pt"
    label = "Simulated post ocean-loss planet orbiting Proxima Centauri"
    sim.smartin.alb_file = "/gscratch/vsm/mwjl/projects/high_res/inputs/desert_highd.alb"
    sim.set_planet_proxima_b()
    sim.load_atmosphere_from_pt(infile, addn2 = False, scaleP = 1.0)
   
    o2 = sim.atmosphere.gases[1]
    o2.cia_file = '/gscratch/vsm/mwjl/projects/high_res/inputs/o4_calc.cia'

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

    sim.smartin.iraylei = [4]
    sim.smartin.vraylei = [1]

    sim.gen_lblscripts()
    sim.run_lblabc()
    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl2 = sim.output.rad.lam
    flux2 = sim.output.rad.pflux
    sflux2 = sim.output.rad.sflux

    adj_flux2 = flux2/sflux2
    return(wl2, flux2)

def ocean_outgassing(lamin, lamax, res):
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

    sim2.lblin.par_file = '/gscratch/vsm/alinc/fixed_input/HITRAN2016.par' #/gscratch/vsm/alinc/fixed_input/
    sim2.lblin.hitran_tag = 'hitran2016'
    sim2.lblin.fundamntl_file = '/gscratch/vsm/alinc/fixed_input/fundamntl2016.dat'
    sim2.lblin.lblabc_exe = '/gscratch/vsm/alinc/exec/lblabc_2016'

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

    sim2.smartin.iraylei = [4]
    sim2.smartin.vraylei = [1]

    sim2.gen_lblscripts()
    sim2.run_lblabc()
    sim2.write_smart(write_file = True)
    sim2.run_smart()

    sim2.open_outputs()
    wl2 = sim2.output.rad.lam
    flux2 = sim2.output.rad.pflux
    sflux2 = sim2.output.rad.sflux

    adj_flux2 = flux2/sflux2 * math.pi
    return(wl2, flux2)

def fluxes(lamin, lamax):
    lamin = 0.76
    lamax = 0.765
    earth_wl, earth_flux = run_earth(lamin,lamax, 0.01)
    wl, flux = run_prox(lamin,lamax, 0.01)
    n_phase = 1000
    phases = np.linspace(0,2*np.pi,n_phase)
    inclination = np.pi/2
    phi_90 = np.pi/2
    sma = 7500000
    c = 299792.458
    
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
    fluxes = []

    while i < n_phase:
        interp = scipy.interpolate.interp1d(obs_wl[:len(flux),i],flux[:len(obs_wl)],fill_value = "extrapolate")
        temp = interp(earth_wl) * earth_flux
        fluxes.append(temp)
        i = i+1
    fluxes = np.asarray(fluxes) 
    fluxes = fluxes.transpose()

    temp = np.arccos(-np.sin(inclination)*np.cos(phases))
    phase_function = ((np.sin(temp)+(np.pi-temp)*(np.cos(temp)))/(np.pi)) 

    fluxes2 = fluxes*(phase_function/phi_90)
    
    from matplotlib.collections import LineCollection

    # Create figure
    fig, ax = plt.subplots(figsize=(12,10))
    ax.set_ylabel("Phase Angle")

    # Create a continuous norm to map from flux to colors
    norm = plt.Normalize(np.min(fluxes2), np.max(fluxes2))

    # Loop over phases
    for i in range(len(phases)):

        # Set dimensions
        x = obs_wl[:,i]
        y = phases[i] * np.ones_like(x)
        z = fluxes2[:,i]

        # Define line points and segments
        points = np.array([x, y]).T.reshape(-1, 1, 2)
        segments = np.concatenate([points[:-1], points[1:]], axis=1)

        # Use linecollections to make color lines
        lc = LineCollection(segments, cmap='bone', norm=norm)
        
        # Set the values used for colormapping
        lc.set_array(z)
        lc.set_linewidth(2)
        line = ax.add_collection(lc)
        
    # Set the axis ranges

    ax.set_xlim(lamin, lamax)
    ax.set_ylim(min(phases), max(phases))

    # Create colorbar
    cbar = fig.colorbar(line)
    cbar.set_label(r"Flux [W/m$^2$/$\mu$m]", rotation = 270, labelpad = 25)



    ax2 = ax.twinx()
    ax2.plot(earth_wl, earth_flux, 'r')
    ax2.set_xlabel(r"Wavelength [$\mu$]")
    ax2.axis('off')
    
    fig.savefig("/gscratch/vsm/mwjl/projects/high_res/plots/" + str(lamin) +  "_RV_thick_line.png", bbox_inches = "tight")

def basic_plot(lamin, lamax, type):
    import platform
    if platform.system() == 'Darwin':
        # On a Mac: usetex ok
        matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})
        matplotlib.rcParams['font.size'] = 25.0
        matplotlib.rc('text', usetex=True)
    elif platform.node().startswith("n"):
        # On hyak: usetex not ok, must change backend to 'agg'
        matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})
        matplotlib.rcParams['font.size'] = 25.0
        matplotlib.rc('text', usetex=False)
        plt.switch_backend('agg')
    earth_wl, earth_flux = run_earth(lamin,lamax, 0.01)
    if type == 0:
         wl, flux = clouds_out(lamin,lamax, 0.01)
    elif type == 1:
         wl, flux = ocean_loss(lamin,lamax, 0.01)
    elif type == 2:
         wl, flux = ocean_outgassing(lamin,lamax, 0.01)
    n_phase = 1000
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
    print(min(rv_bary),max(rv_bary))
    #rv_orb = (np.sqrt((G*m_prox)/(sma_m)) * np.sin(inclination)/1000 *np.sin(phases))
    rv_orb = (sma*2*np.pi)/967680* np.sin(inclination) *np.sin(phases)
    rv_sys = -21.7 * np.ones_like(phases)
    #rv_bary = (29.8 * np.sin((11.2/365.)*phases))
    rv = rv_sys + rv_orb - rv_bary
    inclination = np.pi/2.9
    rv = max(rv)

    obs_wl = np.outer(earth_wl,(1+rv/c))
    fig, ax = plt.subplots(figsize = (10,10))
    ax.plot(earth_wl[:len(earth_flux)], earth_flux[:len(earth_wl)], label = "original")
    ax.plot(obs_wl[:len(flux)], flux[:len(obs_wl)], label = "shifted")
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$m)")
    ax.legend()
    fig.savefig("/gscratch/vsm/mwjl/projects/high_res/plots/" + str(lamin) + str(type) +  "RV_dopp.png", bbox_inches = "tight")

def make_gif(lamin,lamax):
#    if platform.system() == 'Darwin':
#        # On a Mac: usetex ok
#        matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})
#        matplotlib.rcParams['font.size'] = 25.0
#        matplotlib.rc('text', usetex=True)
#    elif platform.node().startswith("n"):
#        # On hyak: usetex not ok, must change backend to 'agg'
#        matplotlib.rc('font',**{'family':'serif','serif':['Computer Modern']})
#        matplotlib.rcParams['font.size'] = 25.0
#        matplotlib.rc('text', usetex=False)
#        plt.switch_backend('agg')
    earth_wl, earth_flux = run_earth(lamin,lamax, 0.01)
    wl, flux = clouds_out(lamin,lamax, 0.01)
    print("wl")
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
    print(min(rv_bary),max(rv_bary))
    #rv_orb = (np.sqrt((G*m_prox)/(sma_m)) * np.sin(inclination)/1000 *np.sin(phases))
    rv_orb = (sma*2*np.pi)/967680* np.sin(inclination) *np.sin(phases)
    rv_sys = -21.7 * np.ones_like(phases)
    #rv_bary = (29.8 * np.sin((11.2/365.)*phases))
    rv = rv_sys + rv_orb - rv_bary
    inclination = np.pi/2.9
    i = 0
    print("loop start")
    while i < len(rv):
        fig, ax = plt.subplots(figsize = (10,10))
        obs_wl = np.outer(temp,(1+rv/c))        
        interp = scipy.interpolate.interp1d(obs_wl[:len(flux),i],flux[:len(obs_wl)],fill_value = "extrapolate")
        temp = interp(earth_wl) * earth_flux
        ax.plot(earth_wl, earth_flux, label = "earth")
        ax.set_ylabel("Reflectance")
        ax.set_xlabel("Wavelength ($\mu$ m)")
        ax.plot(wl, temp, label = "doppler shift")
        ax.legend()
        fig.savefig("/gscratch/vsm/mwjl/projects/high_res/plots/" + str(lamin) + str(i) + "RV.png", bbox_inches = "tight")
        i = i+1

    gif_path = str(lamin) + "RV.gif"
    i = 0
    inputs = []
    while i < len(rv):
        name = "/gscratch/vsm/mwjl/projects/high_res/plots/"+str(lamin)+str(i)+"RV.png"
        inputs.append(name)
        i = i+1
    plt.figure(figsize=(4,4))

    with imageio.get_writer(gif_path, mode='I') as writer:
        for i in range(len(inputs)):
            writer.append_data(imageio.imread(inputs[i].format(i=i)))

def flux_calc(lamin,lamax, type):
    print("flux calc", lamin)
    earth_wl, earth_flux = run_earth(lamin,lamax, high_res)
    earth_wl_low, earth_flux_low = run_earth(lamin, lamax,low_res)
    if type == 0:
         wl, flux = clouds_out(lamin, lamax,high_res)
         wl_low, flux_low = clouds_out(lamin, lamax, low_res)
    elif type == 1:
         wl, flux = ocean_loss(lamin,lamax,high_res)
         wl_low, flux_low = ocean_loss(lamin, lamax, low_res)
    elif type == 2:
         wl, flux = cloud_weight_highw(lamin,lamax, high_res)
         wl_low, flux_low = cloud_weight_highw(lamin, lamax,low_res)

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
    print(len(earth_wl))
    length = np.shape(earth_wl)
    length = length[0]
    out =np.empty((n_phase, length))

    i = 0
    while i < n_phase:
        high_pass_wl, high_pass_flux = high_pass(obs_wl[:len(flux),i],flux[:len(obs_wl)],lamin, lamax,wl_low, flux_low, type) 
        interp = scipy.interpolate.interp1d(high_pass_wl[:len(high_pass_flux)], high_pass_flux[:len(high_pass_wl)],fill_value = "extrapolate")
        temp = np.asarray(interp(earth_wl) * earth_flux)
        out[i,:] = (temp) 
        i = i+1
    return(wl, out)

def high_pass(wl, flux, lamin, lamax, wl_low, flux_low, type):
    long_flux = []
    temp_factor = math.ceil(len(wl) / len(wl_low))
    print(temp_factor)
    for i in flux_low:
        j = 0
        while j < temp_factor: 
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
    return(wl, out)

def integ_calc(lamin,lamax, type):
    print("integ calc", lamin)
    f = open("/gscratch/vsm/mwjl/projects/high_res/scripts/integrations_" + str(lamin) + str(type)  +".txt", "a")
    wl, out = flux_calc(lamin, lamax, type)
    i = 0
    while i <= len(out):
        adds = integrate.trapz(out[:len(wl),i], wl[:len(out)])
        print("adds", lamin, type, adds)    
        name = str(abs(adds)) 
        f = open("/gscratch/vsm/mwjl/projects/high_res/scripts/integrations_" + str(lamin) + str(type)  +".txt", "a")
        f.write(str(name) + "\n")
        f.close()
        i = i+1

def no_phase(lamin, lamax, type):
    if type == 0:
         wl, flux = clouds_out(lamin, lamax,high_res)
         wl_low, flux_low = clouds_out(lamin, lamax, low_res)
    elif type == 1:
         wl, flux = ocean_loss(lamin,lamax, high_res)
         wl_low, flux_low = ocean_loss(lamin, lamax,low_res)
    elif type == 2:
         wl, flux = cloud_weight_highw(lamin,lamax,high_res)
         wl_low, flux_low = cloud_weight_highw(lamin, lamax,low_res)
    wl, out = high_pass(wl, flux, lamin, lamax, wl_low, flux_low, type)
    adds = integrate.trapz(out[:len(wl)], wl[:len(out)])
    print("nophase", lamin, type, adds)
    name = str(abs(adds))
    f = open("/gscratch/vsm/mwjl/projects/high_res/scripts/integrations_nophase.txt", "a")
    f.write(str(name) + "\n")
    f.close()


    
def master_plot_nophase():
    no_phase(0.74, 0.78, 0)
    no_phase(0.67, 0.71, 0)
    no_phase(0.61, 0.65, 0)
    no_phase(1.25, 1.27, 0)
    no_phase(0.74, 0.78, 1)
    no_phase(0.67, 0.71, 1)
    no_phase(0.61, 0.65, 1)
    no_phase(1.25, 1.27, 1)
    no_phase(0.74, 0.78, 2)
    no_phase(0.67, 0.71, 2)
    no_phase(0.61, 0.65, 2)
    no_phase(1.25, 1.27, 2)

def master_plot():
    integ_calc(0.74, 0.78, 0)
    integ_calc(0.67, 0.71, 0)
    integ_calc(0.61, 0.65, 0)
    integ_calc(1.25, 1.27, 0)
    integ_calc(0.74, 0.78, 1)
    integ_calc(0.67, 0.71, 1)
    integ_calc(0.61, 0.65, 1)
    integ_calc(1.25, 1.27, 1)
    integ_calc(0.74, 0.78, 2)
    integ_calc(0.67, 0.71, 2)
    integ_calc(0.61, 0.65, 2)
    integ_calc(1.25, 1.27, 2)

def read_integ():
    print('starting read_integ')
    master_plot_nophase()
    master_plot()

    files = "integrations_0.610.txt","integrations_0.611.txt","integrations_0.612.txt","integrations_0.670.txt","integrations_0.671.txt","integrations_0.672.txt", "integrations_0.740.txt","integrations_0.741.txt","integrations_0.742.txt", "integrations_1.250.txt", "integrations_1.251.txt","integrations_1.252.txt"
    for name in files: 
        print(name)
        output = np.genfromtxt(name)
        print(output)
        phase = np.asarray(range(100))
        phase = phase* 2 * np.math.pi / 100
        print("output, phase", len(output), len(phase))
        phase.astype(np.float)
        output = output[:len(phase)]
        phase = phase[:len(output)]       
    #    fig, ax = plt.subplots(figsize = (12,12))
    #    ax.plot(phase, output[:,1])
    #    ax.set_title("Integration Metric over Phase")
    #    ax.set_ylabel("Integration Metric")
    #    ax.set_xlabel("Phase")
    #    out_name = name[:-4] + ".png"
    #    fig.savefig(out_name, bbox_inches = 'tight')

        integ = integrate.trapz(output,phase) 
        f = open("/gscratch/vsm/mwjl/projects/high_res/scripts/integrations_fin.txt", "a")
        f.write(str(integ) + '\n')
        f.close()
     
if __name__ == '__main__':

    import platform
    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="nor_plt",
                               subname="submit.csh",
                               workdir = "",
                               nodes = 1,
                               mem = "500G",
                               walltime = "108:00:00",
                               ntasks = 28,
                               account = "vsm",
                               submit = True,
                               rm_after_submit = True)
    elif platform.node().startswith("n"):
        # On a mox compute node: ready to run
        print('job submitted') 
#        integ_calc(0.74, 0.78, 0)
#        read_integ()
        flux_calc(0.76, 0.78,0)
    else:
        fluxes(0.60,0.70)

