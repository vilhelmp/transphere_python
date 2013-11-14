









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
        'rin',          20,             'AU',   'float',    # Inner radius of shell
        'rout',         8000,           'AU',   'float',    # Outer radius of shell
        'nshell',       200,            ' ',    'int',      # Number of shells
        'spacing',      'powerlaw1',    ' ',    'str',      # Type of grid spacing
        'nref',         0,              ' ',    'int',      # Refinement shells
        'rref',         0.0,            'AU',   'float',    # Refinement radius
        'rstar',        3.0,            'AU',   'float',    # Stellar radius
        'tstar',        5780,           'K',    'float',    # Stellar temperature
        'mstar',        1,              'MSUN', 'float',    # Stellar mass
        'isrf',         0.0,            ' ',    'mixed',    # Scaling of ISRF
        'tbg',          2.73,           'K',    'float',    # Spectral shape of ISRF (Blackbody equivalent temperature). With -1 the spectrum is read from the file isrf.inp and scaled by ISRF.
        'dpc',          250,            'PC',   'float',    # Distance to source in pc
        'r0',           1.0E3,          'AU',   'float',    # Reference radius
        'plrho',        -1.5,           ' ',    'float',    # Powerlaw for rho
        'rho_type',     'powerlaw1',    ' ',    'str',      # Type of rho dependence with radius ['powerlaw1', 'shu-knee']
        'n0',           2e6,            'cm-3', 'float',    # H2 number density at reference radius
        'r_knee',       1280.0,         'AU',   'float',    # At what radius should the knee be (AU)
        'gas2dust',     100.0,          ' ',    'float',    # Gas to dust ratio
        'localdust',    False,          ' ',    'bool',     # Dust opacity local?
        'silent',       True,           ' ',    'bool',     # Verbose?
        'nriter',       30,             ' ',    'int',      # Maximum nr of iterations
        'convcrit',     1E-5,           ' ',    'float',    # Convergence criterion
        'ncst',         10,             ' ',    'int',      # Nr of rays for star
        'ncex',         30,             ' ',    'int',      # Nr of rays between star and Rin
        'ncnr',         1,              ' ',    'int',      # Nr of rays per radial grid point
        'itypemw',      1,              ' ',    'int',      # Type of mu weighting
        'idump',        True,           ' ',    'bool',     # Dump convergence history?
        'opacfile',     0,              ' ',    'str',      # Filename of opacityfile
        'freqfile',     0,              ' ',    'str',      # Filename of frequency file
        'directory',    '' ,            ' ',    'str'       # Directory to work in
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
        self.r_midpt = array((lower+upper)/2)
        self.r_grid = array([array([i,j]) for i,j in zip(lower,upper)])

        #
        # Create the density array
        #
        # calculate the density at all radial points
        if self.rho_type == 'powerlaw1':
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
            _os.system('rm -f dustopac_1.inp')
            _os.system('rm -f dustopac.inp') # added this, do I really need
                                            # to remove it too?
            f = open(_os.path.join(self.directory, 'dustopac.inp'), 'w')
            f.write('1               Format number of this file'+'\n')
            f.write('1               Nr of dust species'+'\n')
            f.write('============================================================================'+'\n')
            f.write('-1              Way in which this dust species is read (-1=file)'+'\n')
            f.write('0               Extension of name of dustopac_***.inp file'+'\n')
            f.write('----------------------------------------------------------------------------'+'\n')
            f.close
            f = open(_os.path.join(self.directory, 'dustopac_0.inp'), 'w')
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
            _os.system('rm -f {0}'.format(_os.path.join(self.directory, 'dustopac_0.inp')))
            _os.system('rm -f {0}'.format(_os.path.join(self.directory, 'dustopac.inp'))) # added this, do I really need
                                            # to remove it too?
            with open(_os.path.join(self.directory, 'dustopac.inp'),'w') as f:
                f.write('1               Format number of this file\n')
                f.write('1               Nr of dust species\n')
                f.write('============================================================================\n')
                f.write('-1              Way in which this dust species is read (-1=file)\n')
                f.write('1               Extension of name of dustopac_***.inp file\n')
                f.write('----------------------------------------------------------------------------\n')
            with open(_os.path.join(self.directory, 'dustopac_1.inp'), 'w') as f:
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
        with open(_os.path.join(self.directory, 'transphere.inp'),'w') as f:
            f.write(text)
        #
        # Make the stellar information file
        # (mstar and tstar are irrelevant; they are there for historical reasons)
        # isn't "tstar" used for planck calc?
        #~ f=open('starinfo.inp','w')
        with open(_os.path.join(self.directory, 'starinfo.inp'),'w') as f:
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
        with open(_os.path.join(self.directory, 'starspectrum.inp'), 'w') as f:
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

        with open(_os.path.join(self.directory, 'external_meanint.inp'), 'w') as f:
            f.write('{0}\n'.format(len(self.freq)))
            for inu in range(0, len(self.freq)):
                f.write('{0:20}\t{1:20}\n'.format(self.freq[inu], bgspec[inu]))
        #
        # Write the envelope structure
        #
        with open(_os.path.join(self.directory, 'envstruct.inp'),'w') as f:
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
        
