#!/usr/bin/python

####        Modern Earth        ####
# Prepared By:                     #
# Andrew Lincowski                 #
# University of Washington         #
# Virtual Planetary Laboratory     #
####################################

### Python standard libraries ###
import os
from datetime import date
from math import trunc

#### Custom libraries ###
from run_param_multi import *
from build_atm import regrid_atm

# diagonostic 
plot_jacobian = 0 # source 1=solar, 2=thermal, 3=both (have to check function, not fully generalized)

plot_lblabc = 0
lblabc_lyr = 12
lblabc_t = 2
filepath = 'ch4.abs'
wnmin = 700
wnmax = 850

########################
### Run-time options ###
########################

### Run/File handling ###

newrun = 0
runlbl = 0
runvpl = 0
runsmart = 0
runtrans = 0 # transit type: 1) geometric, 2) ray-tracing, 3) MC (warning! takes a long time)
scratch = '/gscratch/vsm/mwjl/scratch'

# Climate
incl_rad = True
incl_conv = True
run_jacobs = 0
toloff = 0
dtemp = [0.5,0.1] # solar, thermal

# Jacobian options for efficiency
split_jacobs = 1 # splits thermal and solar to run in parallel
runsol = 1 # run new solar jacobians
splitsza = 0 # if using multiple SZA's, runs these in parallel (must use split_jacobs and runsol as well)

# Regenerate atm file
genatm = 0
restart = 0
regrid = 0
avgT = 0
smoothT = 0

# Restore previous run
previous = False
prev_date = "_2016-08-13.atm"

# Time-stepping
num_steps = '1000'
num_intervals = '100'
nprint = str(int(num_steps)//200)
dt = '500000'
newflg = True # using the following new methods
itime = 7 # 4 = adaptive pressure-weighted exponential w/inflection, 7 = adaptive, uniform
expp0 = 1e4
dtra = 0.01
dblgrid = True

# Jacobian choices, 0 = none, 1 = TOA rad, 2 = Lyr-by-lyr flux
jacob_pres = 0
jacob_temp = 2 # 2 = flux jacobians (1 = radiance, not for climate use)

# Runs files/folders structure
runfolders = 'earth'
jobstr = 't1e_earth'
smartfile = 'runsmart_' + date.today().isoformat() + '.script'
vplsh = 'runvpl_jt.sh'
outputname = 'vpl_climate'
vplexec = 'vpl_climate_conv79'
if(runsmart):
    vplsh = 'runsmart.sh'
    outputname = 'smart_spectra'
if(run_jacobs == True):
    climatefile = 'runvpl_climate_' + date.today().isoformat() + '.script'
else:
    climatefile = 'runvpl_climate_timestep_' + date.today().isoformat() + '.script'
if(runsmart and not run_jacobs): jacob_temp = 0
restartfile = "restart.atm"

### VPL CLIMATE / SMART BINNING PARAMETERS ###
nstr = '4'
nSZA = '1' #Number of Solar Zenith Angles (computational angles)
tau_err = '0.5' # optical depth error tolerance, lower is better and slower, (0.15 - 0.75)
scat_err = '0.35' # single scattering albedo error tolerance, lower is better and slower, (0.1 - 0.5)
asym_err = '0.35' # scattering asymmetry parameter error tolerance, lower is better and slower, (0.1 - 0.5)
alb_err = '0.25' # surface albedo error tolerance, lower is better and slower, (0.02 - 0.5)
thwnmin = 40.
thwnmax = 5000.
swnmin = 300.
swnmax = 87650.
res = [10.,10.] # climate model fluxes spectral resolution, thermal (1-10), solar (10)
sources = 3 # 1 = solar, 2 = thermal, 3 = both, 4 = transit only

single_gas = 0 # to run multiple SMART instances each with a single gas from the atmosphere
sensitivity = 0 # to run multiple SMART instances each missing one of the atmospheric gases

### SMART Transit parameters ###
n = 1 # yes / no
r_star = 1.
scale_2_solar = 1.
imp_param = 0.

### FIXED MODEL PARAMETERS ###

### PLANET ###
dist = 0.02916 # au
L_yr = 365.25 # planetary days
L_day = 86400. # standard seconds per day
grav = 9.15 # gravity m/s2
radius = 5829. # km

### BULK ATMOSPHERE ###

# Atmospheric structure (P-T; SMART can also use alt-T)
mu_atm = 29.6
atmfile = '../../1bar.atm'
atmskip = 2
atm_cols = [1,2] # Fortran index from 1
nlev = 72

# Rayleigh scattering
atm_codes = [3,4] # 1 = air, 2 = CO2, 3 = N2, 4 = O2, 5 = H2, 6 = He
mix_ratios = [0.79,0.21]

### SURFACE ###
surf_temp = 288.4 # Need to update manually for every run!!
surf_index = 2
albedo_dir = 'albedo/'
surf_files = ['basalt_fresh.alb',
               'blackbody.alb',
               'earth1.alb',
               'earth_desert.alb',
               'earth_m8_turbet.alb',
               'earth_standard.alb',
               'goethite.alb',
               'grey_020.alb',
               'kasting_standard.alb', 
               'mars_avg.alb',
               'melting_snow_15.alb',
               'ocean_ASTER.alb',
               'quartz.alb',
               'snow.alb',
               'snow_fine_ASTER.alb',
               'soil4.alb']
surf_skip = ['16','2','6','6','2','2','102','3','2','3','16','6','15','26','21','26']
surf_cols = '1,2' # Fortran index from 1
surflx = 0. # internal heat flux W/m2

### GASES ###
gas_codes = ([1,7,5,6,3,4,2,22]) # H2O, CO2, O3, N2O, CO, CH4, O2, N2
gas_columns = ([4,5,6,7,8,9,10,11]) # Fortran index from 1
#nfgbrd = ([2,1,2,2,2,2,2,1]) # number of foreign gas broadeners
fb_gas_codes = ([[7,22],[22],[7,22],[7,22],[7,22],[7,22],[7,22],[7]]) # broadener codes
nfgbrd = []
for i in range(len(fb_gas_codes)):
    try:
        nfgbrd.append(len(fb_gas_codes[i]))
    except:
        nfgbrd.append(0.)
num_abs = len(gas_codes)
gas_types = np.array([[1,3],[1,2,3],[1,3],[1,3],[1,3],[1,3],[1,2,3],[1,2,3]]) # types for each gas: 1 = abs, 2 = CIA, 3 = xsec
jacob_gas = [0,0,0,0,0,0,0,0] # 0 = No jacobian, 1 = radiance (not for climate), 2 = flux (for climate)
dmix = [0.,0.,0.,0.,0.,0.,0.,0.] # fractional mixing ratio perturbation
mixing_type = 2 #1 = vol; 2 = mass
mix_scale = [1.,1.,1.,1.,1.,1.,1.,1.]
abs_dir = ''
ntemp = 5 # number of temperature offsets in abs files
tmp_offset = 50. # temperature offsets in abs files

# HITRAN2016
fundfile = 'fundamntl2016.dat'
linefile = 'HITRAN2016.par'
lblexec = '/gscratch/vsm/alinc/exec/lblabc_2016'


### CONDENSIBLES ###
condense = None #[1] # HITRAN gas code(s)
if(condense==None):
    num_cond = 0
else:
    num_cond = len(condense)
conds_skiplines = ([18,23]) # H2O, CO2
conds_cols = ([1,2,3]) # P, T, H; Fortran index from 1
conds_conv = ([100000.,1.]) # convert to Pa, K
reservoir = [1] # Type; 0) dry, 1) fixed humidity, 2) fixed mixing ratio (must then also set if mass or vmr
isat = [2] # 1 = sat trop, 2 = manabe
humid = [70.] # %
readcond = [False]
rstrtskip = 3
rstrt_cols = [2,6,7] #cols of P, mixes # Fortran index from 1

### CONVECTION PARAMETERS ###
conv_type = '2' # [1 = adjustment; 2 = Mixing length scheme; 3 = turbulent; 4 = mix + h2o]
ml_type = '3' # [1 = fixed; 2 = prop to scale height; 3 = Blackadar aymptotic ML]
ml_fix = 1.
ml_prop = 0.01
surf_wind_v = 7.
surf_rough = 0.001

### AEROSOLS ###
jacob_aer = []
dtau = []
aero_codes = None #[]
aero_rad = None #[]
cld_file = fixed_input + '/cld/arbitrary_cld.dat'
aero_wn = None #[]
aero_skip = 3
aero_cols = [] # Fortran index from 1
aero_shgt = []
num_aero = 0# len(aero_rad)
cld_skip = 3
miefile1 = ""
mie_dir = '/gscratch/vsm/alinc/mie_mom/'
miefile = mie_dir + miefile1


# Set aer_str to insert manually-prepared aerosols
# The below are the two standard Earth-model clouds
# Note the Earth model used 50% cloud coverge
# cirrus tau = 5, tau = 20, stratocumulus tau = 5, tau = 10 ?
'''

aer_str = 

2			******NO. AEROSOLS******
0			Aerosol jacobians [0 = None; 1 = Radiance; 2 = Flux]
/gscratch/vsm/alinc/fixed_input/cld/strato_cum
19			Lines to skip, mie
1,7,8,11			Cols of wl, Qext, Qsca, g1
1			Full phase function
17			Lines to skip, mom
/gscratch/vsm/alinc/fixed_input/cld/earth_cld.dat
15400.0			Standard wavenumber [cm^-1]
11			Lines to skip
1			Use pressure coordinate
1,2			columns of P,tau
1.,2.2			Convert to Pa; Scale Tau from 1.0
False			Scale heights?
0			Aerosol jacobians [0 = None; 1 = Radiance; 2 = Flux]
/gscratch/vsm/alinc/fixed_input/cld/baum_cirrus_de100
1			Lines to skip, mie
1,4,5,3			Cols of wl, Qext, Qsca, g1
2			HG
/gscratch/vsm/alinc/fixed_input/cld/earth_cld.dat
15400.0			Standard wavenumber [cm^-1]
3			Lines to skip
1			Use pressure coordinate
1,2			columns of P,tau
1.,3.45			Convert to Pa; Scale Tau from 1.0
False                   Scale heights?'''
aer_str = None

### STAR ###
star = 'TRAPPIST-1'
if(star == 'TRAPPIST-1'):
    if(runsmart):
        star_file = 'specs/TRAPPIST-1_1cm.dat'
    else:
        star_file = 'specs/TRAPPIST-1_10cm.dat'
    star_skip = '10'
    star_cols = '1,2' # Fortran indx from 1
    star_flux_units = '2' #2 = Watts/m**2/micron
    star_spec_units = '1' #1 = Microns; 2 = wn
    star_units2um = '1.0'
elif(star == 'Proxima'):
    if(runsmart):
        star_file = 'specs/proxima_cen.dat'
    else:
        star_file = 'specs/proxima_cen_10.dat'
    star_skip = '25'
    star_cols = '1,2'
    star_flux_units = '2' #2 = Watts/m**2/micron
    star_spec_units = '1' #1 = Microns; 2 = wn
    star_units2um = '1.0'
else: # use the Sun
    star_file = 'specs/Kurucz1cm-1_susim_atlas2.dat'
    star_skip = '12'
    star_cols = '1,2'
    star_flux_units = '1' # W/m2/cm^-1
    star_spec_units = '2' #1 = Microns; 2 = wn
    star_units2um = '1.'

### Deprecated ###
samelbl = False
copyabs = False
multinode = True
runmp = False
wl_bin = 0.001 #um
rundir_hyak = '/scr/mwjl/' + jobstr + '/'

##################################   RUN SCRIPTS   ##################################

orig_dir = os.getcwd()
num_atm = len(atm_codes)

i = 0

# Generate run directory
rundir = runfolders + '/'
if not os.access(orig_dir + '/' + rundir,os.F_OK):
    os.makedirs(orig_dir + '/' + rundir)

# Plotting            
files = [runfolders +  '/' + timefile]
if(plotlast):
    plot_atm(jobstr,files,[plotsteps],nlev,stepper,jobstr + '_' + runfolders+ '_test' + date.today().isoformat(),0) # 1 for h2o
if(plotspec):
    if(split_jacobs or not run_jacobs):
        plot_spec(rundir+outputname+'_timestep',res,[thwnmin,thwnmax],[swnmin,swnmax],nlev,surf_temp,nsza=int(nSZA))
    else:
        plot_spec(rundir+outputname,res,[thwnmin,thwnmax],[swnmin,swnmax],nlev,surf_temp,nsza=int(nSZA))

j = i

os.chdir(rundir)

if(plot_jacobian == 1): # solar
    prefix = 'vpl_climate'
    jsuffix = '_sflx_sza01.j_temp'
    radsuffix = '_toa.r01'
    xran = [0.2,10.0]
    read_jacobian_binary(prefix,jsuffix,radsuffix,swnmin,swnmax,res[1],int(nstr),xran,ncpu=20,pflag='refl')
elif(plot_jacobian == 2): # thermal
    prefix = 'vpl_climate'
    jsuffix = '_tflx.j_temp'
    radsuffix = '_toa.rad'
    xran = [3.0,25.]
    read_jacobian_binary(prefix,jsuffix,radsuffix,thwnmin,thwnmax,res[0],int(nstr),xran,ncpu=20,pflag='flux')
elif(plot_jacobian == 3): # both
    prefix = 'vpl_climate'
    jsuffix = '_sflx_sza01.j_temp'
    radsuffix = '_toa.r01'
    xran = [0.2,3.0]
    read_jacobian_binary(prefix,jsuffix,radsuffix,swnmin,swnmax,res[1],int(nstr),xran,ncpu=20)
    prefix = 'vpl_climate'
    jsuffix = '_tflx.j_temp'
    radsuffix = '_toa.rad'
    xran = [3.0,25.]
    read_jacobian_binary(prefix,jsuffix,radsuffix,thwnmin,thwnmax,res[0],int(nstr),xran,ncpu=20)

if(plot_lblabc):
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    mpl.rcParams.update({'font.size': 14, 'font.family': 'serif'})
    
    P0,temp,wn,absco = read_lblabc(orig_dir + '/' + runfolders+'/' + filepath,wnmin,wnmax)

    print('selected wn',wn[(wn > wnmin) & (wn < wnmax)])
    print('selected abs',absco[lblabc_t][(wn > wnmin) & (wn < wnmax),lblabc_lyr])
    fig,ax = plt.subplots(figsize=[12,6])
    ax.semilogy(wn[(wn > wnmin) & (wn < wnmax)],absco[lblabc_t][(wn > wnmin) & (wn < wnmax),lblabc_lyr])
    plt.savefig('ch4.pdf',bbox_inches='tight')


#G enerate atm files as necessary
if(newrun == True):
     #gen_newrun("clr",orig_dir,rundir,rundir_clr,rundir_AU[i],surf_files[l],atmfile,atmfile0,restartfile)
     os.system('rsync -rptu ' + orig_dir + '/' + atmfile + ' ' + orig_dir + '/' + rundir + atmfile)

if(previous == True): restore_atm(atmfile,restart,previous,prev_date)

if(genatm == True):
    os.chdir(orig_dir)
    gen_atm(restart,rundir,restartfile,atmfile,atmskip,timefile,nlev,gas_columns,atm_cols,mixing_type,mu_atm,gas_codes,usecond=num_cond)

if(regrid):
    p0 = 1. # surface pressure
    nlyrs = 65 # number of layers in new atmosphere # make sure this odd
    pmin3 = 0.01 # top of grid pressure
    pmin2 = 1.0 
    pmin1 = 10.0
    dpmin = 11 # 1 per decade


### Create shell and run scripts ###

if(single_gas): # Single-gas sensitivity tests
    for l in range(num_abs):
        gen_bash_script(orig_dir,rundir,rundir_hyak,vplsh,runvpl,runsmart,runtrans,multinode,runmp,jobstr,runlbl,samelbl,climatefile,smartfile+str(gas_codes[l]),incl_rad,incl_conv,L_day,L_yr,dist,grav,radius,mu_atm,num_atm,atm_codes,atm_names,mix_ratios,jacob_pres,jacob_temp,atmfile,atmskip,atm_cols,ntemp,tmp_offset,[nfgbrd[l]],[fb_gas_codes[l]],surf_temp,num_cond,condense,gas_names,fixed_input,conds_skiplines,conds_cols,conds_conv,1,[gas_codes[l]],[jacob_gas[l]],dmix[l],[gas_types[l]],[gas_columns[l]],mixing_type,[mix_scale[l]],num_aero,jacob_aer,dtau,cld_names,aero_codes,aero_rad,aero_wn,cld_file,cld_skip,aero_cols,aero_shgt,albedo_dir+surf_files[surf_index],surf_skip[surf_index],surf_cols,nstr,star_file,star_cols,star_skip,outputname+str(gas_codes[l]),abs_names,mol_wgt,mie_dir,run_jacobs,num_steps,num_intervals,dt,abs_dir,star_flux_units,star_spec_units,star_units2um,tau_err,scat_err,asym_err,alb_err,dtemp,conv_type,ml_type,ml_fix,ml_prop,surf_wind_v,n,r_star,scale_2_solar,wl_bin,nSZA,vplexec,thwnmin,thwnmax,swnmin,swnmax,copyabs,imp_param,restartfile,rstrtskip,rstrt_cols,readcond,surflx,reservoir,humid,surf_rough=surf_rough,isat=isat,fundfile=fundfile,linefile=linefile,lblexec=lblexec,scratch=scratch,aer_str = aer_str, nprint = nprint)
elif(sensitivity): # Remove-a-gas sensitivity tests
    for l in range(num_abs):
        gen_bash_script(orig_dir,rundir,rundir_hyak,vplsh,runvpl,runsmart,runtrans,multinode,runmp,jobstr,runlbl,samelbl,climatefile,smartfile,incl_rad,incl_conv,L_day,L_yr,dist,grav,radius,mu_atm,num_atm,atm_codes,atm_names,mix_ratios,jacob_pres,jacob_temp,atmfile,atmskip,atm_cols,ntemp,tmp_offset,[x for k,x in enumerate(nfgbrd) if k!=l],[x for k,x in enumerate(fb_gas_codes) if k!=l],surf_temp,num_cond,condense,gas_names,fixed_input,conds_skiplines,conds_cols,conds_conv,num_abs-1,[x for k,x in enumerate(gas_codes) if k!=l],[x for k,x in enumerate(jacob_gas) if k!=l],[x for k,x in enumerate(dmix) if k!=l],[x for k,x in enumerate(gas_types) if k!=l],[x for k,x in enumerate(gas_columns) if k!=l],mixing_type,[x for k,x in enumerate(mix_scale) if k!=l],num_aero,jacob_aer,dtau,cld_names,aero_codes,aero_rad,aero_wn,cld_file,cld_skip,aero_cols,aero_shgt,albedo_dir+surf_files[surf_index],surf_skip[surf_index],surf_cols,nstr,star_file,star_cols,star_skip,outputname,abs_names,mol_wgt,mie_dir,run_jacobs,num_steps,num_intervals,dt,abs_dir,star_flux_units,star_spec_units,star_units2um,tau_err,scat_err,asym_err,alb_err,dtemp,conv_type,ml_type,ml_fix,ml_prop,surf_wind_v,n,r_star,scale_2_solar,wl_bin,nSZA,vplexec,thwnmin,thwnmax,swnmin,swnmax,copyabs,imp_param,restartfile,rstrtskip,rstrt_cols,readcond,surflx,reservoir,humid,parallel=str(gas_codes[l]),ncpu=28,surf_rough=surf_rough,isat=isat,fundfile=fundfile,linefile=linefile,lblexec=lblexec,sources=sources,scratch=scratch,aer_str = aer_str, nprint = nprint)
elif(split_jacobs): # splitting jacobians between solar and thermal
            gen_bash_script_multi(orig_dir,rundir,rundir_hyak,vplsh,runvpl,runsmart,runtrans,multinode,runmp,jobstr,runlbl,samelbl,climatefile,smartfile,incl_rad,incl_conv,L_day,L_yr,dist,grav,radius,mu_atm,num_atm,atm_codes,atm_names,mix_ratios,jacob_pres,jacob_temp,atmfile,atmskip,atm_cols,ntemp,tmp_offset,nfgbrd,fb_gas_codes,surf_temp,num_cond,condense,gas_names,fixed_input,conds_skiplines,conds_cols,conds_conv,num_abs,gas_codes,jacob_gas,dmix,gas_types,gas_columns,mixing_type,mix_scale,num_aero,jacob_aer,dtau,miefile,aero_codes,aero_rad,aero_wn,cld_file,cld_skip,aero_cols,aero_shgt,albedo_dir+surf_files[surf_index],surf_skip[surf_index],surf_cols,nstr,star_file,star_cols,star_skip,outputname,abs_names,mol_wgt,mie_dir,run_jacobs,num_steps,num_intervals,dt,abs_dir,star_flux_units,star_spec_units,star_units2um,tau_err,scat_err,asym_err,alb_err,dtemp,conv_type,ml_type,ml_fix,ml_prop,surf_wind_v,n,r_star,scale_2_solar,wl_bin,nSZA,vplexec,thwnmin,thwnmax,swnmin,swnmax,copyabs,imp_param,restartfile,rstrtskip,rstrt_cols,readcond,surflx,reservoir,humid,toloff=toloff,res=res,runsol=runsol,surf_rough=surf_rough,isat=isat,newflg=newflg,itime=itime,dblgrid=dblgrid,expp0=expp0,dtra=dtra,splitsza=splitsza,alt0=0.,fundfile=fundfile,linefile=linefile,lblexec=lblexec,sources=sources,scratch=scratch,aer_str = aer_str, nprint = nprint)
else: # regular climate / SMART runs
    gen_bash_script(orig_dir,rundir,rundir_hyak,vplsh,runvpl,runsmart,runtrans,multinode,runmp,jobstr,runlbl,samelbl,climatefile,smartfile,incl_rad,incl_conv,L_day,L_yr,dist,grav,radius,mu_atm,num_atm,atm_codes,atm_names,mix_ratios,jacob_pres,jacob_temp,atmfile,atmskip,atm_cols,ntemp,tmp_offset,nfgbrd,fb_gas_codes,surf_temp,num_cond,condense,gas_names,fixed_input,conds_skiplines,conds_cols,conds_conv,num_abs,gas_codes,jacob_gas,dmix,gas_types,gas_columns,mixing_type,mix_scale,num_aero,jacob_aer,dtau,cld_names,aero_codes,aero_rad,aero_wn,cld_file,cld_skip,aero_cols,aero_shgt,albedo_dir+surf_files[surf_index],surf_skip[surf_index],surf_cols,nstr,star_file,star_cols,star_skip,outputname,abs_names,mol_wgt,mie_dir,run_jacobs,num_steps,num_intervals,dt,abs_dir,star_flux_units,star_spec_units,star_units2um,tau_err,scat_err,asym_err,alb_err,dtemp,conv_type,ml_type,ml_fix,ml_prop,surf_wind_v,n,r_star,scale_2_solar,wl_bin,nSZA,vplexec,thwnmin,thwnmax,swnmin,swnmax,copyabs,imp_param,restartfile,rstrtskip,rstrt_cols,readcond,surflx,reservoir,humid,toloff=toloff,res=res,surf_rough=surf_rough,isat=isat,newflg=newflg,itime=itime,dblgrid=dblgrid,expp0=expp0,dtra=dtra,alt0=0.,fundfile=fundfile,linefile=linefile,lblexec=lblexec,sources=sources,scratch=scratch,aer_str = aer_str, nprint = nprint)


