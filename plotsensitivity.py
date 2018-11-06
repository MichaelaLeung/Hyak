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

def longplot(choice):  
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
    if choice == 'h':
        data = smart.readsmart.Rad("longplot/highwsensi_hitran2012_5000_20000cm_toa.rad")
        info = "co2", "o2","o3","co","so2","ocs","n2"
        data1 = smart.readsmart.Rad("longplot/highwsensi_no_o2_hitran2012_5000_20000cm_toa.rad")
        data2 = smart.readsmart.Rad("longplot/highwsensi_no_o3_hitran2012_5000_20000cm_toa.rad")
        data3 = smart.readsmart.Rad("longplot/highwsensi_no_co_hitran2012_5000_20000cm_toa.rad")
        data4 = smart.readsmart.Rad("longplot/highwsensi_no_so2_hitran2012_5000_20000cm_toa.rad")
        data5 = smart.readsmart.Rad("longplot/highwsensi_no_ocs_hitran2012_5000_20000cm_toa.rad")
        data6 = smart.readsmart.Rad("longplot/highwsensi_no_n2_hitran2012_5000_20000cm_toa.rad")
    else:
        data = smart.readsmart.Rad("longplot/highwsensi_hitran2012_11111_11627cm_toa.rad")
        info = "co2", "o2","o3","co","so2","ocs","n2"
        data1 = smart.readsmart.Rad("longplot/highwsensi_no_o2_hitran2012_11111_11627cm_toa.rad")
        data2 = smart.readsmart.Rad("longplot/highwsensi_no_o3_hitran2012_11111_11627cm_toa.rad")
        data3 = smart.readsmart.Rad("longplot/highwsensi_no_co_hitran2012_11111_11627cm_toa.rad")
        data4 = smart.readsmart.Rad("longplot/highwsensi_no_so2_hitran2012_11111_11627cm_toa.rad")
        data5 = smart.readsmart.Rad("longplot/highwsensi_no_ocs_hitran2012_11111_11627cm_toa.rad")
        data6 = smart.readsmart.Rad("longplot/highwsensi_no_n2_hitran2012_11111_11627cm_toa.rad")

    radius = 6850.0
    r_AU = 0.0485

    counter = 0
    for sample in (data, data1, data2, data3, data4, data5, data6):
        wl = sample.lam
        flux = sample.pflux
        sflux = sample.sflux
        adj_flux = flux/sflux * ((radius / r_AU) **2 )
        fig, ax = plt.subplots(figsize = (30, 10))
        ax.plot(wl, adj_flux)
        ax.set_ylabel("Reflectance")
        ax.set_xlabel("Wavelength ($\mu$ m)")
        ax.set_title(info[counter])
        atmos = info[counter]
    fig.savefig("no" + str(atmos) + "highd.png", bbox_inches = 'tight')
    counter = counter+1




if __name__ == '__main__':

    import platform

    if platform.node().startswith("mox"):
        # On the mox login node: submit job
        runfile = __file__
        smart.utils.write_slurm_script_python(runfile,
                               name="sensiplt",
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
        longplot('h')
        counter = 0
    else:
        longplot('m')
        counter = 0








