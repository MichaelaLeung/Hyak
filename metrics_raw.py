import numpy as np
import smart
from matplotlib import pyplot as plt
import matplotlib
from matplotlib.collections import LineCollection
from astropy.io import fits 
import matplotlib as mpl
import sys, os
import datetime
matplotlib.rcParams['text.usetex'] = False
import scipy.integrate as integrate

def run_smart(lamin, lamax):
    HERE = os.path.dirname(os.path.abspath(__file__))
    place = os.path.join(HERE, "metrics")

    try:
        os.mkdir(place)
    except OSError:
        pass

        

    sim = smart.interface.Smart(tag = "prox")
    infile = "profile_Earth_proxb_.pt_filtered"
    res = 1/(15*lamin)
    low_res = 1*lamin
    print(res, low_res)
    sim.set_run_in_place(place) 
    sim.set_executables_automatically()
    sim.set_planet_proxima_b()
    sim.set_star_proxima()
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
    sim.open_outputs()

    wl = sim.output.rad.lam
    flux = sim.output.rad.pflux
    sflux = sim.output.rad.sflux
    adj_flux = flux/sflux * ((sim.smartin.radius / (sim.smartin.r_AU * 149598000)) **2 )
    flux = flux/sflux

    sim2 = smart.interface.Smart(tag = "prox")

    sim2.set_run_in_place(place) 
    sim2.set_executables_automatically()
    sim2.set_planet_proxima_b()
    sim2.set_star_proxima()
    sim2.load_atmosphere_from_pt(infile, addn2 = False)
    sim2.smartin.FWHM = low_res
    sim2.smartin.sample_res = low_res
    sim2.smartin.minwn = 1e4/lamax
    sim2.smartin.maxwn = 1e4/lamin 
    sim2.lblin.minwn = 1e4/lamax
    sim2.lblin.maxwn = 1e4/lamin 
    sim2.gen_lblscripts()
    sim2.run_lblabc()
    sim2.write_smart(write_file = True)
    sim2.run_smart()

 
    sim2.open_outputs()
    wl_low = sim2.output.rad.lam
    flux_low = sim2.output.rad.pflux
    sflux_low = sim2.output.rad.sflux
    flux_low = flux_low/sflux_low


    sim3 = smart.interface.Smart(tag = "earth")

    sim3.set_run_in_place(place) 
    sim3.set_executables_automatically()
    sim3.load_atmosphere_from_pt(infile, addn2 = False)
    sim3.smartin.FWHM = res
    sim3.smartin.sample_res = res
    
    sim3.smartin.minwn = 1e4/lamax
    sim3.smartin.maxwn = 1e4/lamin 
    sim3.lblin.minwn = 1e4/lamax
    sim3.lblin.maxwn = 1e4/lamin 
    sim3.gen_lblscripts()
    sim3.run_lblabc()
    sim3.write_smart(write_file = True)
    sim3.run_smart()

    sim3.open_outputs()
    earth_wl = sim3.output.rad.lam
    earth_flux = sim3.output.rad.pflux
    sflux_earth = sim3.output.rad.sflux
    earth_flux = earth_flux/sflux_earth

    return (wl, flux, adj_flux, wl_low, flux_low, earth_wl, earth_flux)
 
    
def interval(wl, flux):
    lis = []
    st_wl = wl[0]
    delta = []
    exo_dict = dict(zip(flux, wl))
    final_dict = dict()
    cutoff = np.median(flux)
    for i in flux: 
        temp = (i > cutoff)
        if temp == True: 
            start = i 
            st_wl = exo_dict[start]
        elif temp == False: 
            finish = i 
            fin_wl = exo_dict[finish]
            delta = (fin_wl - st_wl)
            final_dict[st_wl] = fin_wl
    return(len(final_dict))

def high_pass(flux, flux_low):
    k = 0
    x = 0
    z = 0
    mixed = []
    out = []
    flattened = []
    long_flux = []
    for i in flux_low:
        j = 0
        while j < 100: 
            long_flux.append(i)
            j = j+1
            
    while k < len(flux):
        temp = (flux[k] + long_flux[k]) / 2
        mixed.append(temp)
        k = k+1
    
    print(len(mixed))
       
    while z < len(mixed): 
        diff = abs(mixed[z] - long_flux[z])
        out.append(diff)
        z = z+1
        
    return(out)



def outputs(lamin, lamax):
    wl, flux, adj_flux, wl_low, flux_low, earth_wl, earth_flux = run_smart(lamin, lamax) #wl, flux, adj_flux, wl_low, flux_low, earth_wl, earth_flux
    adds = max(abs(integrate.trapz((high_pass(flux, flux_low), wl))))
    high = interval(wl, (high_pass(flux, flux_low)))
    label = str(lamin) + "to" + str(lamax)
    out = label, "fpfs", np.median(adj_flux), "line cutoff", high, "integral", adds, "together", (np.median(adj_flux)*high*adds)
    f = open("outputs_small5.txt", "a")
    f.write(str(out) + "\n")

    n_phase = 1000
    phases = np.linspace(0,2*np.pi,n_phase)
    inclination = np.pi/2
    phi_90 = np.pi/2
    sma = 7500000
    c = 299792.458
    
    fluxes = np.outer(adj_flux, np.ones(n_phase)) 
    temp = np.arccos(-np.sin(inclination)*np.cos(phases))
    phase_function = ((np.sin(temp)+(np.pi-temp)*(np.cos(temp)))/(np.pi)) 

    fluxes2 = fluxes*(phase_function/phi_90)

    rv_orb = 35.02* np.sin(inclination) *np.sin(phases)
    rv_sys = -21.7 * np.ones_like(phases)
    rv_bary = (29.8 * np.sin((11.2/365.)*phases))
    rv = rv_sys + rv_orb - rv_bary

    obs_wl = np.outer(wl,(1+rv/c))

    import platform
    if platform.system() == 'Jarvis':
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
    # Create figure
    fig, ax = plt.subplots(2,1, figsize=(20,15))
    ax[0].set_ylabel("Phase Angle")

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
        cmap = mpl.colors.ListedColormap(['cyan',
                                  'royalblue', "rebeccapurple"])
        cmap.set_over('red')
        cmap.set_under('blue')
        bounds = [-0.5, -0.25, 0.0, 0.25, 0.5]
        lc = LineCollection(segments, cmap=cmap, norm=norm)
    
        # Set the values used for colormapping
        lc.set_array(z)
        lc.set_linewidth(2)
        line = ax[0].add_collection(lc)
    
    # Set the axis ranges
    ax[0].set_ylim(min(phases), max(phases))
    ax[0].set_xlim(lamin, lamax)
    

    # Create colorbar
 
    cbar = fig.colorbar(line, pad = 0.05, ax = ax[0], extend = 'both')
    cbar.set_label(r"Reflectance", rotation = 270, labelpad = 25)

    ax[1].plot(wl, adj_flux, 'r')
    ax[1].set_xlabel(r"Wavelength [$\mu$]")

    fig_name = str(lamin) + "to" + str(lamax)
    fig.savefig("plots_raw/" + fig_name +  ".png") 

if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="metrics",
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
        number = range(50,200, 2)
        for i in number:
            i = float(i)
            i = i/100
            outputs(i, i+0.02)
    else:
        # Presumably, on a regular computer: ready to run
       number = range(50,250, 10)
       for i in number:
            i = float(i)
            i = i/100
            outputs(i, i+0.01)
