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

def run_prox(lamin, lamax):

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

def run_earth(lamin, lamax):

    res = 1/(10*lamin)

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
    earth_wl, earth_flux = run_earth(lamin,lamax)
    wl, flux = run_prox(lamin,lamax)
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
    ax.set_xlim(0.76, 0.765)
    ax.set_ylim(min(phases), max(phases))

    # Create colorbar
    cbar = fig.colorbar(line)
    cbar.set_label(r"Flux [W/m$^2$/$\mu$m]", rotation = 270, labelpad = 25)


    fig.savefig("prox0.76a.png", bbox_inches = "tight")
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
    earth_wl, earth_flux = run_earth(lamin,lamax)
    wl, flux = run_prox(lamin,lamax)
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
    ax.plot(wl, flux, label = "original")
    ax.plot(obs_wl, flux, label = "doppler shift")
    ax.set_title(title)
    ax.set_ylabel("Reflectance")
    ax.set_xlabel("Wavelength ($\mu$ m)")
    ax.legend()
    fig.savefig("/gscratch/vsm/mwjl/projects/high_res/plots/" + str(lamin) +  "RV_dopp.png", bbox_inches = "tight")

def make_gif(lamin,lamax)
    mport platform
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
    earth_wl, earth_flux = run_earth(lamin,lamax)
    wl, flux = run_prox(lamin,lamax)
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
    i = 0
    for i < len(rv):
        obs_wl = np.outer(earth_wl,(1+rv/c))
        ax.plot(wl, flux, label = "original")
        ax.plot(obs_wl, flux, label = "doppler shift")
        ax.set_title(title)
        ax.set_ylabel("Reflectance")
        ax.set_xlabel("Wavelength ($\mu$ m)")
        ax.legend()
        fig.savefig("/gscratch/vsm/mwjl/projects/high_res/plots/" + str(lamin) + str(i) + "RV.png", bbox_inches = "tight")
        i = i+1

    gif_path = str(lamin) + "RV.gif"
    inputs = []
    for i< len(rv)
        name = "/gscratch/vsm/mwjl/projects/high_res/plots/"+str(lamin)+str(i)+"RV.png"
        inputs.append(name)
    plt.figure(figsize=(4,4))

    with imageio.get_writer(gif_path, mode='I') as writer:
        for i in range(len(inputs)):
            writer.append_data(imageio.imread(inputs[i].format(i=i)))
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
        fluxes(0.75,0.77)
    else:
        fluxes(0.60,0.70)

