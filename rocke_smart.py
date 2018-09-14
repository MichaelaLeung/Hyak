from __future__ import print_function

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import sys, os
import datetime
import pandas as pd
import multiprocessing
from functools import partial

import smart
from netCDF4 import Dataset

__all__ = ['write_atm', 'SmartRocke']

HERE = os.path.dirname(os.path.abspath(__file__))
FIXED_FILE = os.path.join(os.path.dirname(HERE),
                          "rocke3d_output",
                          "ROCKE3D Gas Profiles_no O3.txt")

MONTHS = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct",
          "nov", "dec"]

def write_atm(atm, columns, file):
        """
        Write an atmospheric structure file (*.atm)

        Parameters
        ----------
        atm : np.ndarray
            2D array of pressure, temperature, and gas mixing ratios
        columns : list
            Names of columns
        file : str
            Name of save file
        """

        formula = columns

        # - fill header string
        header = ''
        for i in range(len(formula)):
            if i == 0:
                header = '{:^9}'.format(str(formula[i]))+(' ')
            elif i == 1:
                header = header + '{:^11}'.format(str(formula[i]))+(' ')
            else:
                header = header + '{:^11}'.format(str(formula[i]))+(' ')

        #write output .atm file to text file
        np.savetxt(file, atm, fmt='%0.5e', header = header, comments=" ")

        return

class SmartRocke(object):

    def __init__(self, ):
        """
        """
        return

    def write_smart_from_rocke3d(self, varfile, statfile, ilon = 0, ilat = 0,
                                 Nlon = 72, Nlat = 46, Npres = 40, tag = "pixel",
                                 write_files = True):
        """

        Parameters
        ----------
        varfile : str
            File name/path for tracked/variable quantities
        statfile : str
            File name/path for untracked/fixed quantities
        ilon : int
            Longitude index
        ilat : int
            Latitude index
        Nlon : int
            Number of longitudes
        Nlat : int
            Number of latitudes
        tag : str
            Prefix name for saved pt file

        Returns
        -------
        ptfile : str
            Name of pt file created
        """

        """
        Read-in and parse tracked variables
        """

        # Specific output file
        data = np.genfromtxt(varfile)

        # Convert mixing ratios ppm --> frac
        data[:,4] = 1e-6*data[:,4]
        data[:,5] = 1e-6*data[:,5]

        columns=["longitude",
                 "latitude",
                 "pressure_mid",
                 "temperature",
                 "O3",
                 "H2O",
                 "liquid water cloud mixing ratio (kg/kg)",
                 "ice cloud mixing ratio (kg/kg)",
                 "ground albedo (%)",
                 "liquid water cloud optical depth",
                 "ice water cloud optical depth"]

        data1 = data.reshape([Nlon, Nlat, Npres, -1])

        """
        Read-in and parse static variables
        """

        # Static atmosphere file
        data = np.genfromtxt(statfile)

        columns=["longitude",
                 "latitude",
                 "pressure_mid"
                 "pressure_top",
                 "N2",
                 "O2",
                 "Ar",
                 "CO2",
                 "CH4"]

        data2 = data.reshape([Nlon, Nlat, Npres, -1])

        """
        Create 2D array for the atmospheric structure
        """

        # Get the pressure and temperature profiles
        Pi = data1[ilon,ilat,:,2]
        Ti = data1[ilon,ilat,:,3]

        # Get the mixing ratios
        Zi = np.vstack([
            data1[ilon,ilat,:,4], # O3
            data1[ilon,ilat,:,5], # H2O
            data2[ilon,ilat,:,4], # N2
            data2[ilon,ilat,:,5], # O2
            data2[ilon,ilat,:,6], # Ar
            data2[ilon,ilat,:,7], # CO2
            data2[ilon,ilat,:,8]  # CH4
                       ])

        # Construct a 2D array with everything
        atm_profiles = np.vstack([Pi, Ti, Zi]).T
        atm_profiles[np.isnan(atm_profiles)] = 1e-12

        # Name columns
        columns = ["Press", "Temp", "O3", "H2O", "N2", "O2", "Ar", "CO2", "CH4"]

        ptfile = tag + ".pt"

        # Write a pt file
        if write_files:
            write_atm(atm_profiles, columns, ptfile)

        """
        Create optical depth files
        """

        # Extract differential optical depths
        icetau = data1[ilon,ilat,:,10]
        wattau = data1[ilon,ilat,:,9]

        # Set nans to zero
        icetau[np.isnan(icetau)] = 0.0
        wattau[np.isnan(wattau)] = 0.0

        # Create 2D array of profiles
        clds = np.vstack([Pi[::-1]*100, wattau[::-1], icetau[::-1]]).T

        # Save file
        cldfile = tag + ".cld"
        columns = ["Press", "Liquid", "Ice"]
        formula = columns
        # - fill header string
        header = ''
        for i in range(len(formula)):
            if i == 0:
                header = '{:^9}'.format(str(formula[i]))+(' ')
            elif i == 1:
                header = header + '{:^11}'.format(str(formula[i]))+(' ')
            else:
                header = header + '{:^11}'.format(str(formula[i]))+(' ')

        # write output
        if write_files:
            np.savetxt(cldfile, clds, fmt='%0.5e', header = header, comments=" ")

        """
        Create surface albedo file
        """

        # Extract value and convert to frac
        alb = 0.01 * data1[ilon,ilat,:,8][0]

        # Wavelength grid for albedo files (SMART will interpolate)
        wl_arr = np.array([0.01, .1, .5, 1.0, 3.0, 3.5])

        # Put albedo value on wavelength grid
        alb_arr = alb * np.ones(len(wl_arr))

        # Set the first and last values to zero (to avoid SMART extrapolation issues)
        alb_arr[0] = 0.0
        alb_arr[-1] = 0.0

        # Create 2D array of arrays
        alb2d = np.vstack([wl_arr, alb_arr]).T

        # Save file
        albfile = tag + ".alb"
        formula = ["wl[um]", "alb"]
        # - fill header string
        header = ''
        for i in range(len(formula)):
            if i == 0:
                header = '{:^9}'.format(str(formula[i]))+(' ')
            elif i == 1:
                header = header + '{:^11}'.format(str(formula[i]))+(' ')
            else:
                header = header + '{:^11}'.format(str(formula[i]))+(' ')

        # write output
        if write_files:
            np.savetxt(albfile, alb2d, fmt='%0.5e', header = header, comments=" ")

        self.ptfile = ptfile
        self.cldfile = cldfile
        self.albfile = albfile
        self.atm2d = atm_profiles
        self.cld2d = clds
        self.alb2d = alb2d

        return

def test_gen_inputs():
    """
    """

    Nlon = 72
    Nlat = 46
    Npres = 40

    # Choose index for longitude and latitude
    ilon = 0
    ilat = 0

    # Choose ROCKE3D output files
    varfile = os.path.join(os.path.dirname(HERE), "rocke3d_output/365day_Rotation/November 365x ROCKE3D Gas Profiles_O3,H2O,clouds.txt")
    fixfile = os.path.join(os.path.dirname(HERE), "rocke3d_output/ROCKE3D Gas Profiles_no O3.txt")

    # Set a place for this run to save files
    place = os.path.join(HERE, "smart_io")

    # Create directory for files, if it doesn't already exist
    try:
        os.mkdir(place)
    except OSError:
        pass

    # Create tag
    tag = os.path.join(place, "pixel")

    inputs = SmartRocke()
    inputs.write_smart_from_rocke3d(varfile, fixfile, ilon = ilon, ilat = ilat,
                                    Nlon = Nlon, Nlat = Nlat, Npres = Npres,
                                    tag = tag, write_files = True)
    ptfile, cldfile, albfile, atm = inputs.ptfile, inputs.cldfile, inputs.albfile, inputs.atm2d

def parse_sims(sim_name = "365day_Rotation", month = "jan"):
    """
    Parse the ROCKE-3D simulations and return the txt and nc files that
    correspond to the simulations you want.

    Parameters
    ----------
    sim_name : str
        Name of the simulation
    month : str
        Abbreviated month to get (e.g. `"jan"`)
    """

    # Path to ROCKE-3D simulations
    simpath = os.path.join(os.path.dirname(HERE), "rocke3d_output")

    # Get list of file and dirs in simulation folder
    ls = os.listdir(simpath)

    unique_sims = []
    selected_sims = []
    for item in ls:
        if "." not in item:
            unique_sims.append(item)
            if sim_name is not None and sim_name.lower() in item.lower():
                selected_sims.append(item)

    if sim_name is None:
        print("`sim_name` options:")
        print(unique_sims)
        return None, None

    if len(selected_sims) > 1:
        print('Multiple simulations match `sim_name`, taking smallest substring')
        lens = np.array([len(item) for item in selected_sims])
        ichoose = np.argmin(lens)
        sim = selected_sims[ichoose]
    elif len(selected_sims) < 1:
        print('No simulations match `sim_name`')
        return None, None
    else:
        sim = selected_sims[0]

    simdir = os.path.join(simpath, sim)

    # Get list of files in simdir
    ls2 = os.listdir(simdir)

    txt_files = []
    nc_files = []
    for item in ls2:
        if item.endswith(".txt"):
            txt_files.append(item)
        if item.endswith(".nc"):
            nc_files.append(item)
        if month.lower() in item.lower():
            if item.endswith(".txt"):
                txt_file = os.path.join(simdir, item)
            if item.endswith(".nc"):
                nc_file = os.path.join(simdir, item)

    return txt_file, nc_file

def run_single_pixel(sim_name, month, ilon, ilat,
                     Nlon = 72, Nlat = 46, Npres = 40):
    """
    Run SMART on a single pixel
    """

    txt_file, nc_file = parse_sims(sim_name = sim_name, month = month)

    # Set a place for this run to save files
    place = os.path.join(HERE, "smart_io2")

    # Create directory for files, if it doesn't already exist
    try:
        os.mkdir(place)
    except OSError:
        pass

    # Create tag
    tag = "%s_%s_pix_%i_%i" %(sim_name, month, ilon, ilat)
    tag_full_path = os.path.join(place, tag)

    # Create smart input files
    inputs = SmartRocke()
    inputs.write_smart_from_rocke3d(txt_file, FIXED_FILE, ilon = ilon, ilat = ilat,
                                    Nlon = Nlon, Nlat = Nlat, Npres = Npres, tag = tag_full_path,
                                    write_files = True)

    ptfile, cldfile, albfile, atm = inputs.ptfile, inputs.cldfile, inputs.albfile, inputs.atm2d

    """
    Do SMART stuff
    """

    # read in pt file to get atm and molwt etc (MUST SCALE PRESSURE CORRECTLY)
    atm = smart.interface.Atmosphere.from_pt(ptfile, atm_dir = place, scaleP=100, addn2=False )

    # Set absolute path to cld template files
    cld_paths = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(smart.__file__))), "fixed_input", "cld")

    # Set specific ice and water cloud mie files
    path_icemie = os.path.join(cld_paths, "baum_cirrus_de100.mie")
    path_watmie = os.path.join(cld_paths, "strato_cum.mie")

    # Define a mie mode for the ice cloud
    mie_ice = smart.interface.MieMode()
    mie_ice.mie_file = path_icemie
    mie_ice.mie_skip = 1           # Lines to skip at top of mie file
    mie_ice.mie_lines = "1,4,5,3"  # Columns of wl, qext, qsca, g1
    mie_ice.iang_smart = 2         # Henyey-Greenestein phase function
    # Define a cloud optical depth for the ice cloud
    cld_ice = smart.interface.CloudTau()
    cld_ice.vert_file = cldfile
    cld_ice.cld_pos = "1,3"
    cld_ice.vert_coord = 1
    cld_ice.vert_skip = 1
    cld_ice.vert_xscale = 1.0
    cld_ice.vert_yscale = 1.0
    cld_ice.vert_ref_wno = 15400.0   # From Earth Model (THIS COULD BE WRONG)
    #cld_ice.vert_ref_wno = 25000.0    # This is 0.4 microns


    # Define a mie mode for the water cloud
    mie_wat = smart.interface.MieMode()
    mie_wat.mie_file = path_watmie
    mie_wat.mie_skip = 17           # Lines to skip at top of mie file
    mie_wat.mie_lines = "1,7,8,11"  # Columns of wl, qext, qsca, g1
    mie_wat.iang_smart = 1          # Full phase function
    # Define a cloud optical depth for the water cloud
    cld_wat = smart.interface.CloudTau()
    cld_wat.vert_file = cldfile
    cld_wat.cld_pos = "1,2"
    cld_wat.vert_coord = 1
    cld_wat.vert_skip = 1
    cld_wat.vert_xscale = 1.0
    cld_wat.vert_yscale = 1.0
    cld_wat.vert_ref_wno = 15400.0   # From Earth Model (THIS COULD BE WRONG)
    #cld_wat.vert_ref_wno = 25000.0    # This is 0.4 microns

    # Combine mie mode and optical depth profiles into an aerosol object
    aerosols = smart.interface.Aerosols(miemodes=[mie_ice, mie_wat],
                                        mietau=[cld_ice, cld_wat])

    # Prep simulation
    sim = smart.interface.Smart(tag = tag, atmosphere = atm, aerosols = aerosols)
    sim.smartin.source = 3     # Solar and thermal
    sim.smartin.out_format = 1 # No transit transmission calculations
    sim.set_run_in_place(place = place)
    sim.set_executables_automatically()

    # Set surface albedo
    sim.smartin.alb_file = albfile
    sim.smartin.alb_skip = 1

    # Run LBLABC
    #sim.lblin.n_t_prof = 5
    #sim.lblin.t_offset = 20.0
    sim.gen_lblscripts()
    sim.run_lblabc()

    # Run SMART
    sim.write_smart()
    sim.run_smart()

    return sim.tag

def func_wrapper(p):
    sim_name, month, ilon, ilat = p[0], p[1], p[2], p[3]
    return run_single_pixel(sim_name, month, ilon, ilat)

def run_test1(processes = 1):
    """
    """

    # Set up multithreading
    pool = multiprocessing.Pool(processes=processes)

    # South pole pixel
    ilon, ilat = 0, 0
    case1 = [["365", mon, ilon, ilat] for mon in MONTHS]

    # Eqatorial Subsolar pixel
    ilon, ilat = 0, 25
    case2 = [["365", mon, ilon, ilat] for mon in MONTHS]

    cases = case1 + case2

    # Run all planets+climates with NIRSpec Prism
    _ = list(
        pool.map(
            func_wrapper,
            cases
        )
    )

    return

def run_test2(processes = 1):
    """
    """

    # Set up multithreading
    pool = multiprocessing.Pool(processes=processes)

    # Create list to append cases to
    cases = []

    # South pole pixel
    ilon, ilat = 0, 0
    cases += [["1d", mon, ilon, ilat] for mon in MONTHS]
    #cases += [["4d", mon, ilon, ilat] for mon in MONTHS]
    #cases += [["8d", mon, ilon, ilat] for mon in MONTHS]
    #cases += [["16d", mon, ilon, ilat] for mon in MONTHS]
    #cases += [["32d", mon, ilon, ilat] for mon in MONTHS]
    cases += [["64d", mon, ilon, ilat] for mon in MONTHS]
    #cases += [["128d", mon, ilon, ilat] for mon in MONTHS]
    #cases += [["256d", mon, ilon, ilat] for mon in MONTHS]
    #cases += [["365d", mon, ilon, ilat] for mon in MONTHS]

    # Eqatorial Subsolar pixel
    ilon, ilat = 0, 22
    cases += [["1d", mon, ilon, ilat] for mon in MONTHS]
    #cases += [["4d", mon, ilon, ilat] for mon in MONTHS]
    #cases += [["8d", mon, ilon, ilat] for mon in MONTHS]
    #cases += [["16d", mon, ilon, ilat] for mon in MONTHS]
    #cases += [["32d", mon, ilon, ilat] for mon in MONTHS]
    cases += [["64d", mon, ilon, ilat] for mon in MONTHS]
    #cases += [["128d", mon, ilon, ilat] for mon in MONTHS]
    #cases += [["256d", mon, ilon, ilat] for mon in MONTHS]
    #cases += [["365d", mon, ilon, ilat] for mon in MONTHS]

    # Run all planets+climates with NIRSpec Prism
    _ = list(
        pool.map(
            func_wrapper,
            cases
        )
    )
    return


if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="smtrock2",
                               subname="submit.csh",
                               workdir = "",
                               nodes = 1,
                               mem = "500G",
                               walltime = "0",
                               ntasks = 28,
                               account = "vsm",
                               submit = True,
                               rm_after_submit = True)
    elif platform.node().startswith("n"):
        # On a mox compute node: ready to run
        run_test2(processes = multiprocessing.cpu_count() - 1)
    else:
        # Presumably, on a regular computer: ready to run
        run_test2(processes = 1)
