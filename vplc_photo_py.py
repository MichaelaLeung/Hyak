#!/usr/bin/python

####        VPLC Photo          ####
# Prepared By:                     #
# Andrew Lincowski                 #
# University of Washington         #
# Virtual Planetary Laboratory     #
####################################

# Reads in.dist, SMART atm files, updates T profile & H2O for atmos

from vplc_photochem import * #import all names from vplc_photochem.py

### Note: in may be necessary to make sure the input PTZ_mixingratios.in file is deleted.
### This code is made to force photochem to use in.dist

### PRIMARY MECHANISMS

genphoto = 0
genvpl = 1

#######################   GENPHOTO   ##########################

# path to atmos in.dist file
atmosfile = './PHOTOCHEM/OUTPUT/out.dist'
#atmosfile = './PHOTOCHEM/in.dist'

# paths to photochem output profile (optional)
ptfile = './PHOTOCHEM/OUTPUT/profile.pt'

# path to VPL/SMART atm files
vplfile = '../1bar.atm'
eddyfile = '../vpl_climate_timestep.conv'
skipatm = 2

# output file for in.dist
outfile = './PHOTOCHEM/in.dist'

ih2o = 2 # column index of H2O in atmos
h2ocol = 3 # column index of H2O in vpl/SMART atm file

# These parameters should match photochem's input .inc file
nz = 200
nq = 71
npar = 2

# Option to Add species
nadd = 0
iadd = None
rmgases =  None #[30,33,44]
nq_out = nq

# Atmosphere construction parameters
# dz should match photgrid.f

mu = 29.76 # molecular weight g/mol
dz = 0.62 # km
newpt = True
alt_max = 124.
alt_min = dz

# eddy diffusion construction parameters
eddy_calc = 0
eddy_read = 1
edd_am = 2.e5
edd_0 = 2.e3
gamma = 0.
grav = 9.8065 # m/s
cp = 1013. # specific heat

# Smooth eddy profile option
# enter # of smoothings to do
smooth_edd = 0

# generate eddy diffusion and T-alt plots
genplots = True

#################################################################

if(genphoto):

    ### Load previous dist file ###
    T,edd,den,o3,denco2,aersol,wfall,rpar,mixin = read_indist(atmosfile,nz,nq,npar)

    # Load VPL Climate profiles and prepare new indist
    mixout,t_out,alt0,nlev,T0,P0 = new_indist(nz,nq,mixin,vplfile,skipatm,PTz_cols=[0,1,2],h2ocol=h2ocol,h2o_mmr=True,mu=mu,ih2o=ih2o,newpt=newpt,ptfile=ptfile,dz=dz,altlbc=0.,nadd=nadd,iadd=iadd,rmgases=rmgases)

    # update eddy diffusion profile
    rho = P0*1e5*mu*1e-3/(8.314*t_out)
    edd = update_eddy(eddyfile,alt0,eddy_read=eddy_read,nlev=nlev,eddy_calc=eddy_calc,smooth_edd=smooth_edd,rho=rho,mu=mu,grav=grav,edd_0=edd_0,edd_am=edd_am)

    ######################
    # manually set values
    #mixout[-2][:] = 1e-38
    #mixout[21][:] = 2e-4
    #####################

    # Write new dist file
    nq = nq_out
    write_indist(t_out,edd,den,o3,denco2,aersol,wfall,rpar,mixout,nz,nq,npar,outfile=outfile)

    # Generate plots
    if(genplots):
        # plot eddy
        fig,ax = plt.subplots(figsize=[8,6])
        ax.semilogx(edd,alt0)
        ax.set_xlabel('Eddy Diffusion [cm^2/s]')
        ax.set_ylabel('Altitude [km]')
        ax.set_ylim(bottom=0.)
        plt.savefig('eddy.pdf',bbox_inches='tight')

        fig,ax = plt.subplots(figsize=[8,6])
        ax.plot(t_out,alt0)
        ax.set_xlabel('Temperature')
        ax.set_ylabel('Altitude [km]')
        plt.savefig('T.pdf',bbox_inches='tight')


#####################################    GENVPL ####################################
'''
#Gas = [output col, HITRAN gas code]
O = [3, 34]
O2 = [4, 7]
H2O = [5, 1]
CO = [11, 5]
CH4 = [18, 6]
C2H6 = [20, 27]
SO2 = [28, 9]
OCS = [34, 19]
O3 = [35, 3]
N2O = [39, 4]
CO2 = [42, 2]
C2H6S = [45, ?] 
HCL = [60, 50]
CH3CL = [62, 24]
'''
ptcols = [3,4,5,11,18,35,39,42]

runtitle = 'testing'
atmfile = '../1bar.atm'
atmskip = 2
N2 = 0.97

photofile = './PHOTOCHEM/OUTPUT/PTZ_mixingratios_out.dist'


# For checking vmix > 1.0
adjnum = False

# Generate new atmfile
genatm = 1
nlyrs = 81
P0 = 1.
nlin = 32

# Move h2o to first col 
lh2o = True
ih2o = 2


if(genvpl):
    photo2vplc(atmfile,atmskip,photofile,ptcols,runtitle,ih2o,N2,genatm,nlyrs,P0,nlin,genplots=0,writeatm=1,lh2o=lh2o)


print('Completed!')
