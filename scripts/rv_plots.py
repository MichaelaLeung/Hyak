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
    return(wl, adj_flux)

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
    return(lam, transmiss)

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

    sim.gen_lblscripts()
    sim.run_lblabc()
    sim.write_smart(write_file = True)
    sim.run_smart()

    sim.open_outputs()
    wl2 = sim.output.rad.lam
    flux2 = sim.output.rad.pflux
    sflux2 = sim.output.rad.sflux

    adj_flux2 = flux2/sflux2
    return(wl2, adj_flux2)

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

def fluxes(lamin, lamax):
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
    wl, flux = run_prox(lamin,lamax, 0.01)
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
    obs_wl = np.outer(earth_wl,(1+rv/c))
    i = 0
    out = []
    while i < len(obs_wl):
        interp = scipy.interpolate.interp1d(obs_wl[:,i],flux,fill_value = "extrapolate")
        out.append(interp(earth_wl) * flux)
        i = i+1
    from matplotlib.collections import LineCollection

    # Create figure
    fig, ax = plt.subplots(figsize=(12,10))
    ax.set_ylabel("Phase Angle")

    # Create a continuous norm to map from flux to colors
    norm = plt.Normalize(np.min(fluxes), np.max(fluxes))
    # Loop over phases
    for i in range(len(phases)):

        # Set dimensions
        x = out[:,i]
        y = phases[i] * np.ones_like(x)
        z = fluxes[:,i]

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
    ax.set_xlim(0.76, 0.765)
    ax.set_ylim(min(phases), max(phases))

    # Create colorbar
    cbar = fig.colorbar(line)
    cbar.set_label(r"Flux [W/m$^2$/$\mu$m]", rotation = 270, labelpad = 25)

    ax2 = ax.twinx()
    ax2.plot(earth_wl, earth_flux, 'r')
    ax2.set_xlabel(r"Wavelength [$\mu$]")

    fig.savefig("/gscratch/vsm/mwjl/projects/high_res/plots/" + str(lamin) +  "_RV.png", bbox_inches = "tight")

def basic_plot(lamin, lamax):
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
    wl, flux = run_prox(lamin,lamax,0.01)
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
    ax.plot(earth_wl, earth_flux, label = "original")
    ax.plot(obs_wl, flux, label = "shifted")
    ax.set_title(title)
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$ m)")
    ax.legend()
    fig.savefig("/gscratch/vsm/mwjl/projects/high_res/plots/" + str(lamin) +  "RV_dopp.png", bbox_inches = "tight")

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
    wl, flux = run_prox(lamin,lamax, 0.01)
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
    earth_wl, earth_flux = run_earth(lamin,lamax, 0.01)
    wl, flux = run_prox(lamin,lamax, 0.01)
    earth_wl_low, earth_flux_low = run_earth(lamin, lamax, 1)
    wl_low, flux_low = run_prox(lamin, lamax, 1)
   
    if type == 0: 
         wl, flux = run_prox(lamin,lamax, 0.01)
         wl_low, flux_low = run_prox(lamin, lamax, 1)
    elif type == 1:
         wl, flux = ocean_loss(lamin,lamax, 0.01)
         wl_low, flux_low = ocean_loss(lamin,lamax, 1)
    elif type == 2:      
         wl, flux = ocean_outgassing(lamin,lamax, 0.01)
         wl_low, flux_low = ocean_outgassing(lamin,lamax, 1)
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
    #rv_orb = (np.sqrt((G*m_prox)/(sma_m)) * np.sin(inclination)/1000 *np.sin(phases))
    rv_orb = (sma*2*np.pi)/967680* np.sin(inclination) *np.sin(phases)
    rv_sys = -21.7 * np.ones_like(phases)
    #rv_bary = (29.8 * np.sin((11.2/365.)*phases))
    rv = rv_sys + rv_orb - rv_bary
    inclination = np.pi/2.9
    obs_wl = np.outer(wl,(1+rv/c))

    fluxes_low = np.outer(earth_flux_low, np.ones(n_phase))
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
    obs_wl_low  = np.outer(wl_low,(1+rv/c))

    i = 0
    while i < n_phase:
        interp = scipy.interpolate.interp1d(obs_wl[:,i],flux,fill_value = "extrapolate")
        out = interp(earth_wl) * earth_flux
        interp_low = scipy.interpolate.interp1d(obs_wl_low[:,i],flux_low,fill_value = "extrapolate")
        out_low = interp(earth_wl_low) * earth_flux_low
        integ_calc(wl, out, wl_low, out_low, i, lamin, type) 
        i = i+1

def integ_calc(wl, flux, wl_low, flux_low, z, lamin, type):
    f = open("/gscratch/vsm/mwjl/projects/high_res/scripts/integrations_" + str(lamin) + "transmiss.txt", "a")
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
    adds = integrate.trapz(out[:len(wl)], wl[:len(out)])
    print(adds)    
    name = str(z) + "   " + str(abs(adds)) 
    f = open("/gscratch/vsm/mwjl/projects/high_res/scripts/integrations_" + str(lamin) + "transmiss" + str(type) +".txt", "a")
    f.write(str(name) + "\n")
    
def master_plot():
    flux_calc(0.74, 0.78, 0)
    flux_calc(0.67, 0.71, 0)
    flux_calc(0.61, 0.65, 0)
    flux_calc(1.25, 1.27, 0)
    flux_calc(0.74, 0.78, 1)
    flux_calc(0.67, 0.71, 1)
    flux_calc(0.61, 0.65, 1)
    flux_calc(1.25, 1.27, 1)
    flux_calc(0.74, 0.78, 2)
    flux_calc(0.67, 0.71, 2)
    flux_calc(0.61, 0.65, 2)
    flux_calc(1.25, 1.27, 2)

def read_integ():
    print('starting read_integ')
#    master_plot()
    files = "integrations_0.61transmiss.txt", "integrations_0.67transmiss.txt", "integrations_0.74transmiss.txt", "integrations_1.25transmiss.txt"
    for name in files: 
        output = np.genfromtxt(name)
        phase = output[:,0] * 2 * np.math.pi / 1000
        phase.astype(np.float)
        fig, ax = plt.subplots(figsize = (12,12))
        ax.plot(phase, output[:,1])
        ax.set_title("Integration Metric over Phase")
        ax.set_ylabel("Integration Metric")
        ax.set_xlabel("Phase")
        out_name = name[:-4] + ".png"
        fig.savefig(out_name, bbox_inches = 'tight')
        integ = integrate.trapz(output[:,1],phase) 
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
                               walltime = "10:00:00",
                               ntasks = 28,
                               account = "vsm",
                               submit = True,
                               rm_after_submit = True)
    elif platform.node().startswith("n"):
        # On a mox compute node: ready to run
        print('job submitted') 
        fluxes(0.75, 0.76)
    else:
        fluxes(0.60,0.70)

