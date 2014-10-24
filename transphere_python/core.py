

import os as os
import sys as sys
import subprocess as subprocess
import scipy as scipy
from scipy import pi
from time import time, sleep
from matplotlib import pyplot as plt; plt.ion()

from .helpers import *

#nice line settings
NiceLineSettings = dict(lw=1, ms=3, mew=0, marker='o')



########################################################################
class Make(object):
    """
    Title
    ------------
    **test**
    Class to interface with Transphere - dust continuum
    radiative transfer code.
    
    
    
    ######################## copied info from old "pyCollapse.py"
    # floats
    rstar    = 5.3453e0 * nc.RS  # Stellar radius
    mstar    = 1 * nc.MS         # Stellar mass
    tstar    = 5000.             # Stellar temperature
    rin      = 200. * nc.AU      # Inner radius of shell
    rout     = 1.2e4 * nc.AU     # Outer radius of shell
    r0       = 1e3 * nc.AU       # Reference radius
    rho0     = 3.46e-18          # Density at reference radius ()
    plrho    = -1.5              # Powerlaw for rho
    
    ### Parameters related to control of code
    
    rref     = 500. * nc.AU      # Refinement radius
    
    # integers
    nr       = 200               # Nr of radial shells for grid
    nref     = 100               # Nr of refinement points
    nriter   = 30                # Maximum nr of iterations
    convcrit = 0.00001           # Convergence criterion
    ncst     = 10                # Nr of rays for star
    ncex     = 30                # Nr of rays between star and Rin
    ncnr     = 1                 # Nr of rays per radial grid point
    itypemw  = 1                 # Type of mu weighting
    
    # booleans
    idump    = 1                 # Dump convergence history
    localdust= 0                 # Dust opacity local?
    
    
    output: 
    
    
    """
    def __init__(self, **kwargs):
        # imports
        from scipy import log10, log, arange, linspace, logspace,\
        array, concatenate, where
        import os
        
        #
        #~ self.silent = silent
        #
        
        
        #~ class Opacity: pass
        """
        localdust : Dust opacity local?
        silent    : Verbose?
        nriter    : Maximum nr of iterations
        convcrit  : Convergence criterion
        ncst      : Nr of rays for star
        ncex      : Nr of rays between star and Rin
        ncnr      : Nr of rays per radial grid point
        itypemw   : Type of mu weighting
        idump     : Dump convergence history
        """
        # Checking input parameters
        params = array([
        'rin',          20,            'AU',   'float',    # Inner radius of shell
        'rout',         8000,          'AU',   'float',    # Outer radius of shell
        'nshell',       200,            ' ',     'int',      # Number of shells
        'spacing',      'powerlaw1',    ' ',     'str',      # Type of grid spacing
        'r_flat',       0,             'AU',   'float',      # Radius where density becomes flat
        't_flat',       0,             'AU',   'float',      # Radius where temperature becomes flat, preceeds t_uplim
        't_uplim',      0,              'K',   'float',    # upper limit of temperature
        'nref',         0,              ' ',     'int',      # Refinement shells
        'rref',         0.0,           'AU',   'float',    # Refinement radius
        'rstar',        3.0,           'AU',   'float',    # Stellar radius
        'tstar',        5780,           'K',   'float',    # Stellar temperature
        'mstar',        1,           'MSUN',   'float',    # Stellar mass
        'isrf',         0.0,            ' ',   'mixed',    # Scaling of ISRF
        'tbg',          2.73,           'K',   'float',    # Spectral shape of ISRF (Blackbody equivalent temperature). With -1 the spectrum is read from the file isrf.inp and scaled by ISRF.
        'dpc',          250,           'PC',   'float',    # Distance to source in pc
        'r0',           1.0E3,         'AU',   'float',    # Reference radius
        'plrho',        -1.5,           ' ',   'float',    # Powerlaw for rho
        'rho_type',     'powerlaw1',    ' ',     'str',      # Type of rho dependence with radius ['powerlaw1', 'shu-knee']
        'n0',           2e6,         'cm-3',   'float',    # H2 number density at reference radius
        'r_knee',       1280.0,        'AU',   'float',    # At what radius should the knee be (AU)
        'gas2dust',     100.0,          ' ',   'float',    # Gas to dust ratio
        'localdust',    False,          ' ',    'bool',     # Dust opacity local?
        'silent',       True,           ' ',    'bool',     # Verbose?
        'nriter',       30,             ' ',     'int',      # Maximum nr of iterations
        'convcrit',     1E-5,           ' ',   'float',    # Convergence criterion
        'ncst',         10,             ' ',     'int',      # Nr of rays for star
        'ncex',         30,             ' ',     'int',      # Nr of rays between star and Rin
        'ncnr',         1,              ' ',     'int',      # Nr of rays per radial grid point
        'itypemw',      1,              ' ',     'int',      # Type of mu weighting
        'idump',        True,           ' ',    'bool',     # Dump convergence history?
        'opacfile',     0,              ' ',     'str',      # Filename of opacityfile
        'freqfile',     0,              ' ',     'str',      # Filename of frequency file
        'directory',    '' ,            ' ',     'str'       # Directory to work in
                    ])
        print ('Model created with the following parameters:')
              
########################################################################
        # need to check that the input is correct.
        # strings should be strings, and floats floats
        class Input: pass
        param_zip = zip(params[0::4],params[1::4], params[2::4], params[3::4])
        for par, stdval, unit, typ in param_zip:
            if par in kwargs: # if input was given, save it
                # print the value that was input for that given parameter
                print '   {0:9} : {1} {2}'.format(par, str(kwargs[par]).strip('\"'), unit)
                value = kwargs[par]
            elif par not in kwargs: # if not input was given, use default
                if stdval != None:
                    # print the default value for the parameter
                    print '   {0:9} : {1:7} {2:4} \t(default)'.format(par, stdval, unit)
                    value = stdval
                else:
                    raise Exception('Wrong parameter input/handling'
                                    'Check input and/or code...')
            # Check what type it should be and store it 
            #
            #string
            if typ == 'str':
                #~ kwargs[par] = '{0}'.format(kwargs[par])
                Input.__dict__[par] = str(value)
                self.__dict__[par] = str(value)
            #integer (not sure about 'itypemw' perhaps boolean)
            elif typ == 'int':
                Input.__dict__[par] = int(value)
                self.__dict__[par] = int(value)
            #float
            elif typ == 'float':
                Input.__dict__[par] = float(value)
                self.__dict__[par] = float(value)
            #bool
            elif typ == 'bool':
                Input.__dict__[par] = bool(value)
                self.__dict__[par] = bool(value)
            # just try to fudge it, see if it works
            elif typ == 'mixed':
                # TODO : change mixed to the same as in the 
                # ratran object
                try:
                    Input.__dict__[par] = int(value)
                    self.__dict__[par] = int(value)
                except ValueError:
                    Input.__dict__[par] = str(value)
                    self.__dict__[par] = str(value)
            else:
                Input.__dict__[par] = value
        
        # input parameters contains all the input needed to 
        # create this class again
        self.Input = Input
        
        # copy important parameters to the main class
        # the parameters in the Input class will never be touched and 
        # kept intact for easy saving/loading
        # perhaps create it as a hidden class, i.e. start with "_"?
        #~ self.rin = self.Input.rin
        #~ self.rout = self.Input.rout
        #~ self.rref = self.Input.rref
        #~ self.nref = self.Input.nref
        #~ self.nshell = self.Input.nshell
        #~ self.gas2dust = self.Input.gas2dust
        #~ self.n0 = self.Input.n0        
        #~ self.r0 = self.Input.r0   
        #~ self.rstar = self.Input.rstar
        #~ self.mstar = self.Input.mstar
        #~ self.isrf = self.Input.isrf
        #~ self.tbg = self.Input.tbg
        #~ self.tstar = self.Input.tstar
        #~ self.dpc = self.Input.dpc
        
        #~ self.directory = self.Input.directory
        
        # CHECK if directory exists 
        # create dir if it doesn't exist
        # will raise error if things go south 
        # (i.e., permissions not correct [I think...])
        input_dir_path = os.path.join(os.getcwd(), self.directory)
        if not make_dirs(input_dir_path): 
            print('Directory exists, continuing.')
        
        # copy the frequency file to the directory where we are 
        # running transphere in the end
        _subprocess.Popen(['cp', self.freqfile, self.directory])
        
        # save the transphere input
        # TODO : move outside of the class so its optional?
        save_transphere(self)
        
        # Convert distances to CGS units
        self.rin   *= _cgs.AU   # cm
        self.rout  *= _cgs.AU   # cm
        self.rref  *= _cgs.AU   # cm
        self.r0    *= _cgs.AU   # cm
        self.rstar *= _cgs.RSUN # cm
        self.mstar *= _cgs.MSUN # g
        self.r_flat *= _cgs.AU  # cm
        self.t_flat *= _cgs.AU  # cm
        
        if self.rref:
            # if we want refinement, just call the grid function twice
            # using the self.rin, self.rref and self.nref for first grid
            # and self.rref and self.rout and self.nshell for second
            from scipy import concatenate
            print('Creating grid with refinement.')
            inner_grid = create_grid(
                        self.rin, self.rref, self.nref,
                        space = self.spacing, end=True
                                    )
            outer_grid = create_grid(
                        self.rref, self.rout, self.nshell,
                        space = self.spacing, end=True
                                    )
            radii = concatenate([inner_grid[:-1], outer_grid])
            #~ return inner_grid, outer_grid, radii
        else:
            # if we dont want refinement
            radii = create_grid(self.rin, 
                                self.rout, 
                                self.nshell,
                                space = self.spacing, 
                                end = True)
        #
        # radii is now in cm
        # separate the different radii, the lower/upper bounds
        # and the mid-point of the grid "cells"
        lower = radii[:-1]
        upper = radii[1:]
        self.radii = radii
        self.r_grid = array([array([i,j]) for i,j in zip(lower,upper)])

        #
        # Create the density array
        #
        # calculate the density at all radial points
        if self.rho_type == 'powerlaw1':
            
            # this is dependent on the 
            self.r_midpt = array((lower+upper)/2)
            # 1E-2  - gas-to-dust
            # get the DUST density
            #~ self.rho_gas = self.rho0_gas * (self.radii / self.r0)**(self.plrho)
            r_dependence = (self.radii / self.r0)**(self.plrho)
            # now also save the dust density, and H2 number density
            #~ self.rho_gas = self.rho0_gas * (self.radii / self.r0)**(self.plrho)
            #~ self.n_h2 = self.rho_gas / (_cgs.MUH2 * _cgs.MP)
            #~ self.rho_dust = self.n_h2 * 1 / self.gas2dust
            #
            # n0 is in nh2 / cm3
            #    g/cm3           nh2 / cm3 * amu/h2 * mass of amu(g)
            #~ self.rho0_gas = self.n0 * _cgs.MUH2 * _cgs.MP   # g * cm-3
            # all gas, MUH2 is including H2+He+Metals
            # apply the radial dependence to the gas density
            #~ self.rho_gas = self.rho0_gas * r_dependence     # g * cm-3
            # apply radial dependence to the H2 number density
            self.n_h2 = self.n0 * r_dependence              # nh2 * cm-3
            self.rho_gas = self.n_h2 * _cgs.MUH2 * _cgs.MP  # g * cm-3
            self.rho_dust = self.rho_gas / self.gas2dust  # g * cm-3
            # gas to dust is mass relationship
            # g * cm-3 
            #~ self.rho_dust = self.rho0_gas / self.Input.gas2dust * r_dependence 

        elif self.rho_type == 'shu_knee':
            
            self.r_midpt = array((lower+upper)/2)
            # where to put the knee
            self.r_knee = self.r_knee * _cgs.AU # in cm
            #
            # n(r) = n0 (r / r0) **-1.5     for r <= r_knee
            # n(r) = n0 (r / r0) **-2.0     for r >= r_knee
            #
            #
            # two r_dependence
            # get the n0 for each r_dependence
            self.n01 = self.n0
            self.n02 = self.n0 * (self.r_knee / self.r0)**-1.5 / (self.r_knee / self.r0)**-2.0
            
            # get the inner/outer radii cells
            cells_in = where(self.radii <= self.r_knee)
            cells_out = where(self.radii > self.r_knee)
            r_dep_inner = (self.radii[cells_in] / self.r0)**(-1.5)
            r_dep_outer = (self.radii[cells_out] / self.r0)**(-2.0)
            
            # self.n_h2 = self.n0 * r_dependence
            n_inner = self.n01 * r_dep_inner
            n_outer = self.n02 * r_dep_outer
            # put them in the same array
            self.n_h2 = concatenate((n_inner, n_outer))                 # cm-3
            # calculate the gas density
            self.rho_gas = self.n_h2 * _cgs.MUH2 * _cgs.MP              # g * cm-3
            self.rho_dust = self.rho_gas /self.gas2dust                 # g * cm-3

            # apply radial dependence to the H2 number density
                          # nh2 * cm-3
            # gas to dust is mass relationship
            # g * cm-3 
            #~ self.rho_dust = self.rho0_gas / self.Input.gas2dust * r_dependence 
            
            #r_dependence = concatenate((r_dep_inner, r_dep_outer))
            
            #
            # now also save the dust density, and H2 number density
            #~ self.rho_gas = self.rho0_gas * (self.radii / self.r0)**(self.plrho)
            #~ self.n_h2 = self.rho_gas / (_cgs.MUH2 * _cgs.MP)
            #~ self.rho_dust = self.n_h2 * 1 / self.gas2dust
            #
            # n0 is in nh2 / cm3
            #    g/cm3           nh2 / cm3 * amu/h2 * mass of amu(g)
            #~ self.rho0_gas = self.n0 * _cgs.MUH2 * _cgs.MP   # g * cm-3
            # all gas, MUH2 is including H2+He+Metals
            # apply the radial dependence to the gas density
            #~ self.rho_gas = self.rho0_gas * r_dependence     # g * cm-3
            # apply radial dependence to the H2 number density
                          # nh2 * cm-3
            # gas to dust is mass relationship
            # g * cm-3 
            #~ self.rho_dust = self.rho0_gas / self.Input.gas2dust * r_dependence 
        elif self.rho_type == 'flat_powerlaw1':
            # 1E-2  - gas-to-dust
            # get the DUST density
            #~ self.rho_gas = self.rho0_gas * (self.radii / self.r0)**(self.plrho)
            r_dependence = (self.radii / self.r0)**(self.plrho)
            ind = (self.radii <= self.r_flat).nonzero()[0]
            if len( ind ) < 1:
                raise Exception('No points inside of r_flat')
            else:
                r_dependence[ind] = r_dependence[ind[-1]]
            # now also save the dust density, and H2 number density
            #~ self.rho_gas = self.rho0_gas * (self.radii / self.r0)**(self.plrho)
            #~ self.n_h2 = self.rho_gas / (_cgs.MUH2 * _cgs.MP)
            #~ self.rho_dust = self.n_h2 * 1 / self.gas2dust
            #
            # n0 is in nh2 / cm3
            #    g/cm3           nh2 / cm3 * amu/h2 * mass of amu(g)
            #~ self.rho0_gas = self.n0 * _cgs.MUH2 * _cgs.MP   # g * cm-3
            # all gas, MUH2 is including H2+He+Metals
            # apply the radial dependence to the gas density
            #~ self.rho_gas = self.rho0_gas * r_dependence     # g * cm-3
            # apply radial dependence to the H2 number density
            self.n_h2 = self.n0 * r_dependence              # nh2 * cm-3
            self.rho_gas = self.n_h2 * _cgs.MUH2 * _cgs.MP  # g * cm-3
            self.rho_dust = self.rho_gas / self.gas2dust  # g * cm-3
            # gas to dust is mass relationship
            # g * cm-3 
            #~ self.rho_dust = self.rho0_gas / self.Input.gas2dust * r_dependence 
        else:
            raise Exception('rho_type parameter not recognised')

        
        ################################################################        
        from scipy import array
        from sys import exit as sysexit
        #~ import os
        # if a dust opacity file is given
        # that is it is not tabulated correctly
        # as a dustopac.inp file
        if self.opacfile not in [0, '0']:
            if self.freqfile in [0, '0']:
                sysexit('No frequency file given (freqfile). Cannot proceed. '
                'Need two files: one with frequency, and one with opacity'
                ' at corresponding frequency.')
            with open(self.opacfile, 'r') as f:
                lines = f.read().split('\n')
            # get the "header" info,
            # that is the number of tabulated entries in
            # absorption and scattering each.
            try:
                nf, ns = [int(i) for i in lines[:2][0].split()]
            except:
                errmsg = 'Error parsing Opacity-file. Wrong format.'
                raise Exception(errmsg)
            if nf >= 1:
                # shape of the arrays, for absorption and scattering
                # (nf, ns)
                # so we have to nest the list comprehensions, and
                # call float() on the smallest possible element
                # (i.e. when ns > 1).
                try:
                    cabs = array([[float(i) for i in j.split()] for j in lines[2 : nf+2]])
                    csca = array([[float(i) for i in j.split()] for j in lines[nf+2 : 2*nf+2]])
                except:
                    errmsg = 'Error parsing Opacity-file. Wrong format.'
                    raise Exception(errmsg)
                nrtrange = 1
                trange = [0.e0,0.e0]
            else:
                ########################################################
                # this part not edited, dont have a file like this
                sysexit('This part of the code is not tested'
                        'Please optimize code and check it.')
                nf = ns
                ismooth = 0
                nrtrange = 0
                with open(self.opacfile, 'r') as f:
                    line = f.readline().strip()
                    ns, ismooth, nrtrange = line.split()
                    ns, ismooth, nrtrange = int(ns), int(ismooth), int(nrtrange)
                    if ismooth != 0:
                        import sys
                        sys.exit('Error: Smoothing not yet allowed.')
                    import numpy as np
                    cabs = np.zeros((nf, ns), float)
                    csca = np.zeros((nf, ns, nrtrange), float)
                    dum = 0.e0
                    trange = np.zeros(nrtrange+1, float)
                    
                    for ir in range(0, nrtrange):
                        a, b = f.readline().strip().split()
                        a, b = int(a), int(b)
                        trange[ir] = b
                        for kk in range(0, nf):
                            for iss in range(0, ns):
                                dum = float(f.readline().split())
                                cabs[kk, iss, ir] = dum
                        for kk in range(0,nf):
                            for iss in range(0, ns):
                                dum = float(f.readline().split())
                                csca[kk, iss, ir] = dum
                ########################################################
            with open(self.freqfile, 'r') as f:
                lines = f.read().split('\n')
            nf = int(lines[0])
            freq = array(lines[2:2+nf], 'float')
            wave = 2.9979e14 / freq

            ## TODO : change this, do I even need a dictionary?
            # move the calls to where the variable is defined
            class Opacity: pass
            self.Opacity = Opacity
            
            self.Opacity.ns         = ns
            self.Opacity.nf         = nf
            self.Opacity.freq       = freq
            self.Opacity.wave       = wave
            self.Opacity.cabs       = cabs
            self.Opacity.csca       = csca
            self.Opacity.nrtrange   = nrtrange
            self.Opacity.trange     = trange
            
            self.freq               = freq
            #~ self.opacity = {'ns': ns, 'nf': nf, 'freq': freq, 'wave': wave, 'cabs': cabs, 'csca': csca, 'nrt': nrtrange, 'trange': trange}
        else:
            # always neeed a opacity file and a frequency file
            sysexit('No opacity-file given.')
        ################################################################
        # if local dust opacity
        # the opacity is tabulated with radius
        #
        # have not got a file like this, so i am not changing it

        if self.localdust:
            print ('variable localdust not False/0, this part of the code'
            'is not up to date.')
            sysexit('This part of the code is not tested'
                    'Please optimize code and check it.')
            if 'nr' not in kwargs:
                sysexit('Number of radial shells, \"nr\" not given as input.')
            os.system('rm -f dustopac_1.inp')
            os.system('rm -f dustopac.inp') # added this, do I really need
                                            # to remove it too?
            f = open(os.path.join(self.directory, 'dustopac.inp'), 'w')
            f.write('1               Format number of this file'+'\n')
            f.write('1               Nr of dust species'+'\n')
            f.write('============================================================================'+'\n')
            f.write('-1              Way in which this dust species is read (-1=file)'+'\n')
            f.write('0               Extension of name of dustopac_***.inp file'+'\n')
            f.write('----------------------------------------------------------------------------'+'\n')
            f.close
            f = open(os.path.join(self.directory, 'dustopac_0.inp'), 'w')
            f.write(nr)
            f.write(str(nf)+' 1\n')
            f.write(' ')
            redux=1.e0
            
            for ir in range(0,nr):
                for inu in range(0,nr):
                    f.write(self.Opacity.cabs[inu] * redux)
                for inu in range(0,opacity['nf']):
                    f.write(self.Opacity.csca[inu] * redux)
                f.write(' ')
            f.close
        elif not self.localdust:
            # first remove the standard ratran dust opacity input files
            os.system('rm -f {0}'.format(os.path.join(self.directory, 'dustopac_0.inp')))
            os.system('rm -f {0}'.format(os.path.join(self.directory, 'dustopac.inp'))) # added this, do I really need
                                            # to remove it too?
            with open(os.path.join(self.directory, 'dustopac.inp'),'w') as f:
                f.write('1               Format number of this file\n')
                f.write('1               Nr of dust species\n')
                f.write('============================================================================\n')
                f.write('-1              Way in which this dust species is read (-1=file)\n')
                f.write('1               Extension of name of dustopac_***.inp file\n')
                f.write('----------------------------------------------------------------------------\n')
            with open(os.path.join(self.directory, 'dustopac_1.inp'), 'w') as f:
                f.write(str(self.Opacity.nf)+' 1\n \n')
                for inu in range(0, self.Opacity.nf):
                    f.write(str(self.Opacity.cabs[inu][0])+'\n')
                for inu in range(0, nf):
                    f.write(str(self.Opacity.csca[inu][0])+'\n')
    
    def write_transphereinput(self):
        from scipy import pi, zeros
        import os
        # Transphere input file
        text = ('{0}\n{1}\n{2}\n{3}\n'
                '{4}\n{5}\n{6}\n{7}'.format(2,
                self.nriter,
                self.convcrit,
                self.ncst,
                self.ncex,
                self.ncnr,
                self.itypemw,
                int(self.idump)))
        with open(os.path.join(self.directory, 'transphere.inp'),'w') as f:
            f.write(text)
        #
        # Make the stellar information file
        # (mstar and tstar are irrelevant; they are there for historical reasons)
        # isn't "tstar" used for planck calc?
        #~ f=open('starinfo.inp','w')
        with open(os.path.join(self.directory, 'starinfo.inp'),'w') as f:
            f.write('1\n'
                    '{0}\n'
                    '{1}\n'
                    '{2}\n'.format(self.rstar,
                                self.mstar,
                                self.tstar))
        #
        # The stellar spectrum
        #
        #~ f=open('starspectrum.inp','w')
        sspec = (self.rstar / _cgs.PC)**2 * pi * bplanck(self.freq, self.tstar)
        with open(os.path.join(self.directory, 'starspectrum.inp'), 'w') as f:
            f.write('{0}\n'.format(len(self.freq)))
            for inu in range(0,len(self.freq)):
                f.write('{0:20}\t{1:20}\n'.format(self.freq[inu], sspec[inu]))
        #
        # The exterior spectrum
        #
        if self.tbg == 0.0 or self.isrf == 0:
            bgspec = zeros((len(self.freq)),float)
        elif self.tbg == -1:
            f = open('isrf.inp', 'r')
            nf = int(f.readline().strip())
            bgspec = zeros((len(self.freq)),float)
            for ii in range(0,nf):
                bgspec[ii] = float(f.readline().strip())*self.isrf
            f.close()
        else:
            if self.tbg > 0: bgspec = bplanck(self.freq, self.tbg) * self.isrf

        with open(os.path.join(self.directory, 'external_meanint.inp'), 'w') as f:
            f.write('{0}\n'.format(len(self.freq)))
            for inu in range(0, len(self.freq)):
                f.write('{0:20}\t{1:20}\n'.format(self.freq[inu], bgspec[inu]))
        #
        # Write the envelope structure
        #
        with open(os.path.join(self.directory, 'envstruct.inp'),'w') as f:
            f.write(str(len(self.radii))+'\n')
            f.write(' '+'\n')
            # rho_dust has unit g/cm3
            for ir in range(0,len(self.radii)):
                f.write("%13.6E %13.6E %13.6E" % (self.radii[ir], self.rho_dust[ir], 0.e0)+'\n') # ,format='(3(E13.6,1X))'
    
    def run(self, nice = 0):
        import subprocess
        from time import time, sleep
        import sys
        
        #
        # re-check the input files here?
        # at least the format?
        #
        print ('Running Transphere')
        # temporarily change directory, if it has been initiated
        # the directory should exist
        with ChangeDirectory(self.directory): 
            t1 = time()
            if nice:
                proc = subprocess.Popen(['nice', 'transphere'], 
                                    stdout = subprocess.PIPE, 
                                    stderr = subprocess.STDOUT)
            elif not nice:
                proc = subprocess.Popen(['transphere'], 
                                    stdout = subprocess.PIPE, 
                                    stderr = subprocess.STDOUT)
            sys.stdout.write('Iteration no : ')
            sys.stdout.flush() # flush output so we can write again
            trans_out = []
            while True:
                # first : if process is done, break the loop
                if proc.poll() != None: 
                    break
                nextline = proc.stdout.readline()
                trans_out.append(nextline)
                
                if "Iteration" in nextline:
                    # grab the iteration number
                    iter_no = int(nextline[10:])
                    sys.stdout.write('{0:3} \b\b\b\b'.format(iter_no)) # print out a point
                    sys.stdout.flush()
                    
                #~ if "Error" in nextline:
                    # if it is the error of the first iteration
                    # we grab the error of it
                    #~ if iter_no == 1:
                        #~ start_error = float(nextline[13:])
                        #~ first = 0 # now we are not at the first any longer
                    # if it is not the first error, calculate per cent done
                    #~ else:
                        #~ current_error = float(nextline[13:])
                        #~ diff1 = start_error - self.convcrit
                        #~ diff2 = start_error - current_error
                        #~ p_done = (diff2 / diff1 * 100)
                        #~ sys.stdout.write(' {0} : {1}\r'.format(iter_no, p_done)) # print out a point
                        #~ sys.stdout.flush()    # flush output so we can write again
                #~ sleep(0.5)            # wait for 0.5 second
                
                #~ sys.stdout.write('Downloading File FooFile.txt [%d%%]\r'%i)
                #~ ys.stdout.flush()

            print('\nDone in {0:2.1f} seconds'.format((time()-t1)))
        # now get the output as a string, so we can check stuff if
        # we want to
        self.transphere_output = ''.join(trans_out)
        # read in the output-files, for checks, plots etc
        self.Envstruct, self.Convhist = read_transphereoutput(self)
        if self.t_uplim:
            ind = (self.Envstruct.temp >= self.t_uplim).nonzero()[0]
            self.Envstruct.temp[ind] = self.t_uplim
        elif self.t_flat:
            ind = (self.radii <= self.t_flat).nonzero()[0]
            #~ ind = (self.Envstruct.temp >= self.t_uplim).nonzero()[0]
            self.Envstruct.temp[ind] = self.Envstruct.temp[max(ind)]
        
class Transphere(object):
    
    
    def __init__(self, modelfile = 'trans_input.inp'):
        # read in the model parameters and convert some of them
        # to relevant types
        model, settings = read_transphere_modelfile(modelfile)
        # create the wavelength grid
        freq, lam = create_wavelength_grid(
                    wl_interval = model['wl_interval'], 
                    wl_points =  model['wl_points'])
        # read the opacity table
        raw_opacities = get_opac(infile = model['opacfile'])
        # interpolate opacity table to frequency 
        interp_opacities = interpolate_opacities(
                opac = raw_opacities, 
                freq = freq)
        # find kappa at 550 nm (why?)        
        kappa = find_kappa(raw_opacities['lam'], raw_opacities['kabs'], lam0=550.0, localdust = model['localdust'])
        # create radial grid        
        r = transphere_make_r_grid(
            rin = model['rin'], 
            rout = model['rout'], 
            rref = model['rref'], 
            nout = model['nout'], 
            nin = model['nin'])
        # calculate density profile
        n_h2, rho_gas, rho_dust = transphere_calculate_density(
            r, 
            n0 = model['n0'], 
            r0= model['r0'], 
            plrho = model['plrho'], 
            g2d = model['g2d'],  # default is 100, is set in reader
            rho0 = None)
        # calculate tau
        # kappa is dust opacities, so rho_dust should be input, right?
        tau = transphere_calculate_tau(r, rho_dust, kappa)
        # stellar spectrum
        #~ sspec=(model['rstar'] /nc.pc)**2*math.pi*astroProcs.bplanck(model['freq'],model['tstar'])
        stellar_spec = (float(model['rstar']) * co.R_sun.cgs.value / co.pc.cgs.value)**2. * pi * plancknu(freq, float(model['tstar']))*1e-23     
        # store some things as attributes
        self.model = model
        self.settings = settings
        self.freq = freq 
        self.lam = lam
        self.raw_opacities = raw_opacities
        self.interp_opacities = interp_opacities
        self.kappa = kappa
        self.r = r 
        self.n_h2 = n_h2
        self.rho_gas = rho_gas
        self.rho_dust = rho_dust
        self.tau = tau
        self.stellar_spec = stellar_spec
        self.directory = settings['dirname']
        return None
        
    def write_input(self):
        """
        
        Writes all input to a sub-directory provided in the 
        setup file.
        
        """
        # create dir
        input_dir_path = os.path.join(os.getcwd(), self.directory)
        if not make_dirs(input_dir_path): 
            print('Directory exists, continuing.')

        # Dustopacities        
        with open(os.path.join(self.directory, 'dustopac.inp'),'w') as f:
            f.write('1               Format number of this file\n')
            f.write('1               Nr of dust species\n')
            f.write('============================================================================\n')
            f.write('-1              Way in which this dust species is read (-1=file)\n')
            f.write('1               Extension of name of dustopac_***.inp file\n')
            f.write('----------------------------------------------------------------------------\n')
            f.close()
        
        # write the dustopac_1 file
        #TODO: need to update the writing here...
        with open(os.path.join(self.directory, 'dustopac_1.inp'),'w') as f:
            f.write('{0}  1 \n\n'.format(len(self.freq)))
            for inu in range(len(self.freq)):
                f.write('%13.6f\n'%(self.interp_opacities['kabs'][inu]))
            for inu in range(len(self.freq)):
                f.write('%13.6f\n'%(0.0))
            f.close()
    
        # Transphere settings input file        
        text = ('{0}\n{1}\n{2}\n{3}\n'
                '{4}\n{5}\n{6}\n{7}'.format(2,
                self.settings['nriter'],
                self.settings['convcrit'],
                self.settings['ncst'],
                self.settings['ncex'],
                self.settings['ncnr'],
                self.settings['itypemw'],
                int(self.settings['idump'])))
        with open(os.path.join(self.directory, 'transphere.inp'),'w') as f:
            f.write(text)
        
        # Frequency input file
        text = ['{0:13.7e}\n'.format(f) for f in self.freq ]
        with open(os.path.join(self.directory, 'frequency.inp'),'w') as f:
            f.write('{0} \n\n'.format(len(self.freq)))
            f.writelines(text)

        # Central star information 
        with open(os.path.join(self.directory, 'starinfo.inp'),'w') as f:
            # CGS units!
            f.write('1\n'
                '{0:13.7e}\n'
                '{1:13.7e}\n'
                '{2:13.8f}\n'.format(float(self.model['rstar']) * co.R_sun.cgs.value,
                            float(self.model['mstar']) * co.M_sun.cgs.value,
                            float(self.model['tstar']) ))
    
        # Stellar spectrum in the wavelength range, black body emission
        # 'plancknu' gives output in Jy, multiply with 1e-23 cgs units
        # erg / (s cm^2 Hz)
        text = ['{0:13.7e}\t{1:13.7e}\n'.format(i,j) for i,j in zip(self.freq, self.stellar_spec)]
        with open(os.path.join(self.directory, 'starspectrum.inp'), 'w') as f:
            f.write('{0} \n'.format(len(self.freq)))
            f.writelines(text)
        
        # interstellar radiation field
        if float(self.model['tbg']) == 0.0 and not self.model.has_key('isrf'):
            bgspec = zeros((len(self.freq)),float)
        elif self.model.has_key('isrf'):
            # read the isrf file in and store in bgspec
            f = open(self.model['isrf'], 'r')
            nf = int(f.readline().strip())
            bgspec = zeros((len(self.freq)),float)
            for ii in range(0,nf):
                bgspec[ii] = float(f.readline().strip())
            f.close()
        elif self.model['tbg'] > 0: 
                bgspec = plancknu(self.freq, float(self.model['tbg']) ) * 1e-23
        text = ['{0:13.7e}\t{1:13.7e}\n'.format(i,j) for i,j in zip(self.freq, bgspec)]
        with open(os.path.join(self.directory, 'external_meanint.inp'), 'w') as f:
            f.write('{0}\n'.format(len(self.freq)))
            f.writelines(text)

        with open(os.path.join(self.directory, 'envstruct.inp'),'w') as f:
            f.write('{0} \n\n'.format(str(len(self.r))))
            # rho_dust has unit g/cm3
            for ir in range(0,len(self.r)):
                f.write("%13.6E %13.6E %13.6E\n" % (self.r[ir] * co.au.cgs.value, self.rho_dust[ir], 0.e0))
       
        return None

    def run(self, nice = 0):

        
        #
        # re-check the input files here?
        # at least the format?
        #
        print ('Running Transphere')
        # temporarily change directory, if it has been initiated
        # the directory should exist
        with ChangeDirectory(self.directory): 
            t1 = time()
            if nice:
                proc = subprocess.Popen(['nice', 'transphere'], 
                                    stdout = subprocess.PIPE, 
                                    stderr = subprocess.STDOUT)
            elif not nice:
                proc = subprocess.Popen(['transphere'], 
                                    stdout = subprocess.PIPE, 
                                    stderr = subprocess.STDOUT)
            sys.stdout.write('Iteration no : ')
            sys.stdout.flush() # flush output so we can write again
            trans_out = []
            while True:
                # first : if process is done, break the loop
                if proc.poll() != None: 
                    break
                nextline = proc.stdout.readline()
                trans_out.append(nextline)
                
                if "Iteration" in nextline:
                    # grab the iteration number
                    iter_no = int(nextline[10:])
                    sys.stdout.write('{0:3} \b\b\b\b'.format(iter_no)) # print out a point
                    sys.stdout.flush()
                    
                #~ if "Error" in nextline:
                    # if it is the error of the first iteration
                    # we grab the error of it
                    #~ if iter_no == 1:
                        #~ start_error = float(nextline[13:])
                        #~ first = 0 # now we are not at the first any longer
                    # if it is not the first error, calculate per cent done
                    #~ else:
                        #~ current_error = float(nextline[13:])
                        #~ diff1 = start_error - self.convcrit
                        #~ diff2 = start_error - current_error
                        #~ p_done = (diff2 / diff1 * 100)
                        #~ sys.stdout.write(' {0} : {1}\r'.format(iter_no, p_done)) # print out a point
                        #~ sys.stdout.flush()    # flush output so we can write again
                #~ sleep(0.5)            # wait for 0.5 second
                
                #~ sys.stdout.write('Downloading File FooFile.txt [%d%%]\r'%i)
                #~ ys.stdout.flush()

            print('\nDone in {0:2.1f} seconds'.format((time()-t1)))
        # now get the output as a string, so we can check stuff if
        # we want to
        self.transphere_output = ''.join(trans_out)
        # read in the output-files, for checks, plots etc
        #~ self.Envstruct, self.Convhist = read_transphereoutput(self)
        #~ if self.t_uplim:
            #~ ind = (self.Envstruct.temp >= self.t_uplim).nonzero()[0]
            #~ self.Envstruct.temp[ind] = self.t_uplim
        #~ elif self.t_flat:
            #~ ind = (self.radii <= self.t_flat).nonzero()[0]
            #ind = (self.Envstruct.temp >= self.t_uplim).nonzero()[0]
            #~ self.Envstruct.temp[ind] = self.Envstruct.temp[max(ind)]

# run transphere    except (ZeroDivisionError):
# read in results        return 0
# plot results

def read_transphereoutput(self, ext = 0):
    from scipy import array
    if ext == 0: ext=''
    filetoread = 'envstruct' + ext + '.dat'
    path_to_file = os.path.join(self.directory, filetoread)
    with open(path_to_file, 'r') as f:
        lines = f.read().split('\n')
    nr = int(lines[0])
    dat_envstruct = array([i.split() for i in lines[2:nr + 2]], dtype='float')


    filetoread = 'spectrum.dat'
    path_to_file = os.path.join(self.directory, filetoread)
    with open(path_to_file, 'r') as f:
        lines = f.read().split('\n')
    nr = int(lines[0])
    dat_spectrum = array([i.split() for i in lines[2:nr+2]], dtype='float')

    #~ nr = int(f.readline().strip())
    #~ dat = np.zeros((3,nr),float)
    #~ for ii in range(0,3):
        #~ for jj in range(0,nr):
            #~ dum = f.readline().strip()
            #~ dat[ii,jj] = dum
    #~ f.close()
    class Envstruct:
        r = dat_envstruct[:,0]
        rho_dust = dat_envstruct[:,1]
        temp = dat_envstruct[:,2]

    #~ class Spectrum:
    Envstruct.frequency = dat_spectrum[:,0]
    Envstruct.intensity = dat_spectrum[:,1]
    #~ Envstruct.Spectrum = Spectrum
    #~ self.Envstruct = Envstruct

    #~ import numpy as np
    filetoread = 'convhist.info'
    path_to_file = os.path.join(self.directory, filetoread)
    f = open(path_to_file, 'r')
    nn = int(f.readline().strip().split()[0])
    f.close()

    # Convergence history
    filetoread = 'convhist.dat'
    path_to_file = os.path.join(self.directory, filetoread)
    with open(path_to_file, 'r') as f:
        lines = f.read().split('\n')
    nr = int(lines[0].strip())
    if nr == 0: raise Exception('Nothing run, no convergence history.')
    x1 = nr+1

    #These need to depend on value of nr
    dat1 = array([i.split() for i in lines[1:x1]], dtype='float')
    dat2 = array([i.split() for i in lines[x1+1:x1*2]], dtype='float')
    dat3 = array([i.split() for i in lines[x1*2+1:x1*3]], dtype='float')

    dat = array([dat1,dat2,dat3])

    #~ f = open('convhist.dat','r')
    #~ nr = int(f.readline().strip())
    #~ dat = np.zeros((9,nn,nr),float)
    #~ for jj in range(0,nn):
        #~ for kk in range(0,nr):
            #~ dum = f.readline().strip().split()
            #~ if dum == []: dum=f.readline().strip().split()
            #~ dat[0:9,jj,kk]=np.array(dum,dtype=float)
    #~ f.close()

#    if nn gt 1 then idx=[1,2,0] else idx=[1,0]. Note transpose commands not executed...
    class Convhist:
        temp=dat[:,:,0]
        jjme=dat[:,:,1]
        hhme=dat[:,:,2]
        jj=  dat[:,:,3]
        hh=  dat[:,:,4]
        kapt=dat[:,:,5]
        kapj=dat[:,:,6]
        kaph=dat[:,:,7]
        fj=  dat[:,:,8]
    #~ self.Convhist = Convhist

    #~ f = open('envstruct.inp')
    #~ nr = int(f.readline().strip())
    #~ dat = np.zeros((3,nr),float)
    #~ for ii in range(0,nr):
        #~ dum=f.readline().strip().split()
        #~ if dum == []: dum=f.readline().strip().split()
        #~ dat[0:3,ii]=np.array(dum,dtype=float)
    #~ r=dat[0,:]
    #~ f.close()

    #~ convhist={'r': r, 'temp': temp, 'jjme': jjme, 'hhme': hhme, 'jj': jj, 'hh': hh, 'kapt': kapt, 'kapj': kapj, 'kaph': kaph, 'fj': fj}
    #~ self.Envstruct = envstruct
    #~ self.convhist = convhist
    if self.t_uplim:
        ind = (Envstruct.temp >= self.t_uplim).nonzero()[0]
        Envstruct.temp[ind] = self.t_uplim
    elif self.t_flat:
        ind = (self.radii <= self.t_flat).nonzero()[0]
        #~ ind = (self.Envstruct.temp >= self.t_uplim).nonzero()[0]
        Envstruct.temp[ind] = Envstruct.temp[max(ind)]
    return Envstruct, Convhist

def create_grid(r_in, r_out, nshell, space = 'powerlaw1', end = True):
    # function to create grid
    if space == 'log10':
        from scipy import log10, logspace
        # get the exponent of the start- and
        # stop-radius in input units
        start = [log10(r_in), 0][r_in == 0]
        stop = log10(r_out)
        radii = logspace(start, stop, num=nshell, endpoint=end)
    elif space in ["powerlaw1", "flat_powerlaw1"]:
        from scipy import arange
        radii = r_in * (r_out/r_in)**(arange(nshell)/(nshell - 1.0))
    elif space == 'linear':
        from scipy import linspace
        # linearly spaced grid
        radii = linspace(r_in, r_out, num=nshell, endpoint=end)
    elif space == 'powerlaw2':
        from scipy import linspace
        # first check if coefficients to the power-law was given
        #~ if 'exp' in kwargs:
            #~ p_exp = kwargs['exp']
        #~ else: # if not, set it to 2, i.e. r^2
            #~ p_exp = 2
        radii = r_in + (r_out - r_in)*(linspace(r_in, r_out, num=nshell, endpoint=end)/(r_out))**2
        #pr_int('Not implemented yet.')
        #raise ParError(spaced)
    else:
        raise Exception(space)
    return radii

# FIXME, does not work tries to import old adavis module
def plot_spectrum(freq, intensity, dpc = 0, jy = 0, pstyle = '', xlog = 1, ylog = 1):
    import sys
    import matplotlib.pyplot as pl
    pl.ion()
    #~ from ..views import set_rc
    #~ set_rc
    xcoord = 1.0e4 * _cgs.CC / freq

    if dpc == 0: sys.exit('Error: distance needs to be set when plotting flux')

    distfact = 1.e0/ (dpc**2)

    if jy != 0:
        lumfact = 1e+23
    else:
        lumfact = freq

    pl.plot(xcoord, distfact * lumfact * intensity, pstyle)
    pl.xlabel(r'$\lambda\, [\mu \mathrm{m}]$')

    if jy != 0:
        pl.ylabel(r'$F_\nu$\, [Jy]')
    else:
        pl.ylabel(r'$\nu F_\nu \, [\mathrm{erg cm}^{-2}\, \mathrm{s}^{-1}]$')

    if xlog == 1: pl.xscale('log')
    if ylog == 1: pl.yscale('log')


# FIXME, does not work tries to import old adavis module

# temporary function
# needs to be more modular
def plot_envstruct(self, mol_abundance = '', mark100k = True, **kawargs):
    if not hasattr(self, 'Envstruct'):
        raise Exception('you havent read in the transphere output')
    import matplotlib.pyplot as pl
    from matplotlib.ticker import ScalarFormatter, LogFormatter
    pl.ion()
    pl.close()
    fig = pl.figure(1, figsize=(8,6))
    ax1 = fig.add_subplot(111)
    pl.grid()
    ax2 = ax1.twinx()
    # Density
    p1 = ax1.loglog(self.Envstruct.r/_cgs.AU, self.n_h2, label='n_H2', **kawargs)
    ax1.set_xlabel('Radius (AU)')
    #~ ax1.set_xscale('log')
    ax1.set_ylabel('Number Density (cm-3)')
    # Temperature
    p2 = ax2.loglog(self.Envstruct.r/_cgs.AU, self.Envstruct.temp, color='r', label='Temp', **NiceLineSettings)

    ax2.yaxis.set_major_formatter(ScalarFormatter())
    
    ax2.set_ylabel('Temp (K)', color='r')
    if mol_abundance != '':
        def make_patch_spines_invisible(ax):
            ax.set_frame_on(True)
            ax.patch.set_visible(False)
            for sp in ax.spines.itervalues():
                sp.set_visible(False)
        ax3 = ax1.twinx()
        #~ ylims = ax3.get_ylim()
        #ax3.set_ylim(-0.05E-7, 1.85E-7)

        #~ p3 = ax3.loglog(self.Envstruct.r/_cgs.AU, mol_abundance, 'g')
        p3 = ax3.semilogx(self.Envstruct.r/_cgs.AU, mol_abundance, color='g', label='Mol Abund', **NiceLineSettings)
        ax3.spines["right"].set_position(("axes", 1.2))
        make_patch_spines_invisible(ax3)
        ax3.spines["right"].set_visible(True)
        ax3.set_ylabel('Rel. Abundance', color='g')
        #ax1.legend([p1, p2, p3], ['Density', 'Temp', 'Rel. Abundance'])
        #~ ax3.yaxis.set_major_formatter(())
        #~ ax3.xticks(['1E-9','1E-8','1E-7','1E-6','1E-5','1E-4'],[1E-9,1E-8,1E-7,1E-6,1E-5,1E-4])
        #~ ax3.set_yticks([1E-9,1E-8,1E-7,1E-6,1E-5,1E-4], minor=True)
        #~ ax3.tick_params(axis='y', direction='in')
        fig.subplots_adjust(right = 0.75)
    
    if mark100k:
        from scipy import where
        # where is the value closest to 100 K?
        i_100k = where(abs(100 - self.Envstruct.temp).round(2) == round(min(abs(100 - self.Envstruct.temp)), 2))[0][0]
        r_100k = self.Envstruct.r[i_100k]/_cgs.AU
        t_100k = self.Envstruct.temp[i_100k]
        ax2.annotate('T = {0:.1f} K\nR = {1:.1f} AU'.format(round(t_100k,2), r_100k),
                xy=(r_100k, t_100k), xycoords='data',
                xytext=(-30, -100), textcoords='offset points', fontsize=12,
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=-.2"))
        #~ ax2.plot([r_100k], [t_100k] , 'o',color='r', ms=4, mew=0)
        #~ pl.legend('n_H2', 'Temp', 'Mol Abund')
    #~ else:
    ax1.xaxis.set_major_formatter(ScalarFormatter())
    if mol_abundance == '':
        #Create custom artists
        simArtist = pl.Line2D((0,1),(0,0), color='b')
        anyArtist = pl.Line2D((0,1),(0,0), color='r')
        
        #Create legend from custom artist/label lists
        ax1.legend([simArtist,anyArtist],
                  ['Density', 'Temperature'])
    elif mol_abundance != '':
        #Create custom artists
        simArtist = pl.Line2D((0,1),(0,0), color='b')
        anyArtist = pl.Line2D((0,1),(0,0), color='r')
        molArtist = pl.Line2D((0,1),(0,0), color='g')
        
        #Create legend from custom artist/label lists
        ax1.legend([simArtist, anyArtist, molArtist],
                  ['Density', 'Temperature', 'Mol. abundance'])
    
def cleanup_transphere():
    import os
    filelist = ['convhist.info', 'external_meanint.inp' , 'spectrum.dat', 'transphere.dat', 'dustopac_1.inp',  'envstruct.dat', 'starinfo.inp', 'transphere.inp', 'convhist.dat',  'dustopac.inp', 'envstruct.inp', 'starspectrum.inp']
    for f in filelist:
        os.system('rm {0}'.format(f))

def save_transphere(Obj, filename = 'transpheremodel.pickle'):
    # take input object
    # and save it as a dictionary to filename in directory
    import pickle
    # take all the original input parameters and put into a dictionary
    inp = vars(Obj.Input)
    # now, a dictionary is a non dynamic structure
    # as opposed to a dynamically created object attribute
    # e.g., with the __dict__ method (-> doesn't work with pickle)
    with ChangeDirectory(Obj.directory):
        with open(filename, 'w') as f:
            pickle.dump(inp, f)

def load_transphere(directory = '', filename = 'transpheremodel.pickle'):
    # load input object from filename in directory
    # and create the transphere object
    import pickle
    with open(os.path.join(directory, filename), 'r') as f:
        inputdict = pickle.load(f)
    Obj = Transphere(**inputdict)
    # IDEA : add so that it loads the output(?) as well?
    return Obj

def _pdfcheck(pdf):
        if not pdf: 
            _plt.ion()
        elif pdf:
            _plt.ioff()

def _pdfsave(pdf, pdfname, **kwargs):
    if pdf:
            _plt.savefig('{0}.pdf'.format(str(pdfname)), kwargs)

