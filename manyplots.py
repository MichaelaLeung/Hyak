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

def manyplot(res, lamin, lamax):

    sim = smart.interface.Smart(tag = "prox")
    sim2 = smart.interface.Smart(tag = "earth")
    infile = "profile_Earth_proxb_.pt_filtered"
    earth_infile = "earth_avg.pt"
    
    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, "manyplot")

    try:
        os.mkdir(place)
    except OSError:
        pass

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

    sim2.set_run_in_place(place) 
    sim2.set_executables_automatically()
    sim2.load_atmosphere_from_pt(infile, addn2 = False)
    sim2.smartin.FWHM = res
    sim2.smartin.sample_res = res
    print(lamax)
    
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

    sim2.open_outputs()
    earth_wl = sim2.output.rad.lam
    earth_flux = sim2.output.rad.pflux

    n_phase = 1000
    phases = np.linspace(0,2*np.pi,n_phase)
    inclination = np.pi/2
    phi_90 = np.pi/2
    sma = 7500000
    c = 299792.458
 
    
    flux2 = []
    for i in flux: 
        if i > np.median(flux): 
            flux2.append(np.median(flux))
        else: 
            flux2.append(i)

    fluxes = np.outer(flux2, np.ones(n_phase)) 
    temp = np.arccos(-np.sin(inclination)*np.cos(phases))
    phase_function = ((np.sin(temp)+(np.pi-temp)*(np.cos(temp)))/(np.pi)) 

    fluxes2 = fluxes*(phase_function/phi_90)

    rv_orb = 35.02* np.sin(inclination) *np.sin(phases)
    rv_sys = -21.7 * np.ones_like(phases)
    rv_bary = (29.8 * np.sin((11.2/365.)*phases))
    rv = rv_sys + rv_orb - rv_bary

    obs_wl = np.outer(wl,(1+rv/c))
    import platform
    if platform.system() == 'Darwin':
        # On a Mac: usetex ok
        mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})
        mpl.rcParams['font.size'] = 25.0
        mpl.rc('text', usetex=True)
    elif platform.node().startswith("n"):
        # On hyak: usetex not ok, must change backend to 'agg'
        mpl.rc('font',**{'family':'serif','serif':['Computer Modern']})
        mpl.rcParams['font.size'] = 25.0
        mpl.rc('text', usetex=False)
        plt.switch_backend('agg')
    
    scaled_flux = []
    for i in earth_flux: 
        if i > (max(earth_flux)/5): 
            scaled_flux.append(max(earth_flux))
        else: 
            scaled_flux.append(i)

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
    ax.set_ylim(min(phases), max(phases))

    # Create colorbar
    cbar = fig.colorbar(line)
    cbar.set_label(r"Flux [W/m$^2$/$\mu$m]", rotation = 270, labelpad = 25)

    ax2 = ax.twinx()
    ax2.plot(earth_wl, earth_flux, 'r')
    ax2.set_xlabel(r"Wavelength [$\mu$]")

    fig_name = str(lamin) + "to" + str(lamax)
    fig.savefig(fig_name +  ".png")
 

if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="manyplt",
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
        number = range(50,200,1)
        for i in number:
            i = float(i)
            i = i/100
            manyplot(0.01, i, i+0.1)
    else:
        # Presumably, on a regular computer: ready to run
        manyplot(1, 0.5,0.51)








