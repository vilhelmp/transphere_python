

import os as _os
#~ import sys as _sys
#~ import subprocess as _subprocess
#~ import scipy as _scipy
#~ from matplotlib import pyplot as _plt





########################################################################
# GENERAL HELP FUNCTIONS (move to adavis_core)
def check_input(input_dictionary, input_defaults):
    return 0

def make_dirs(path):
    import errno
    # the makedirs function will raise a EEXIST error 
    # if the directory already exists, if so it returns False
    # if some other error is raise, it will raise that 
    try:
        _os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            return False
    return True


class ChangeDirectory:
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = _os.getcwd()
        _os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        _os.chdir(self.savedPath)

# Now you can enter the directory like this:

#~ with cd("~/Library"):
    #~ # we are in ~/Library
    #~ run some code
    #~ import subprocess
    #~ subprocess.call("ls")





import scipy as sp
import numpy as np
import adapy.libs.cgsconst as _cgs
import matplotlib.pyplot as plt; plt.ion()
import os as os
import astropy.units as un
import astropy.constants as co

OPACDIR = "opacities"

pc2au = 206264.806 # 1 PC in AU, or arcsecond in ???
# 206264.806 is also the number of arcsecs in a circle / 2pi

deg2rad = lambda deg : deg * sp.pi / 180

rad2deg  = lambda rad : rad * 180 / sp.pi

def minmax(arr):
    return [arr.min(), arr.max()]

def plancknu(nu, temp):
    """
    Calculate the planck function at a frequency and temperature
    Temperature can be an array.
    
    Input
    -----
    nu : frequency (Hz)
    temp : temperature (Kelvin) 
    
    Output
    ------
    Jansky (Jy)
    
    """
    import adapy.libs.cgsconst as _cgs
    import scipy as _sp
    import numpy as _np
    x =(_cgs.HH*nu)/(_cgs.KK*temp)
    #bbnu = 2.0*hh*nu*nu*nu/(cc*cc)*(1.0/(np.exp( (hh*nu)/(kk*temp)) - 1.0))
    try:    # assume temp is an array
        bbnu = sp.ones_like(temp)
        bbnu = 2.0 * _cgs.HH * nu**3/(_cgs.CC**2)*(1.0/(sp.exp(x) - 1.0))
        # make sure we correct the function the same as if INT
        bbnu[x == 0] = 0.0
        ind = x < 1e-5 
        bbnu[ind] = 2.0 * nu**2*_cgs.KK*temp[ind]/(_cgs.CC**2)
        ind = x > 20.0
        bbnu[ind] = 2.0 * _cgs.HH * nu**3/(_cgs.CC**2) * sp.exp(-1.0*x[ind])
    except (TypeError): # if TypeError is raised, its an integer
                        # so treat it like one
        if x == 0:
            bbnu = 0.0
        elif x < 1e-5:
            bbnu = 2.0 * nu**2*_cgs.KK*temp/(_cgs.CC**2)
        elif x > 20.0:
            bbnu = 2.0 * _cgs.HH * nu**3/(_cgs.CC**2) * sp.exp(-1.0*x)
        else:
            bbnu = 2.0 * _cgs.HH * nu**3/(_cgs.CC**2)*(1.0/(sp.exp(x) - 1.0))
    # return Jy
    return bbnu*1e23

def plancknu_units(nu, temp):
    """
    Calculate the planck function at a frequency and temperature
    Temperature can be an array.
    
    Input
    -----
    nu : frequency (Hz)
    temp : temperature (Kelvin) 
    
    Output
    ------
    Astropy Quantity of Jansky (Jy)
    
    
    """
    #~ import adapy.libs.cgsconst as _cgs
    # new imports
    from scipy import exp, ones_like
    import astropy.constants as co
    import astropy.units as un
    #~ print nu, temp
    try:
       tmp1 = nu.unit 
       tmp2 = temp.unit
    except (AttributeError):
        print ('Need input units of frequency and temperature.')
        print ('Please use astropy.units module to provide units')
        return None
    #~ x =(_cgs.HH*nu)/(_cgs.KK*temp)
    x =(co.h.cgs*nu)/(co.k_B.cgs*temp) # dimensionless
    #bbnu = 2.0*hh*nu*nu*nu/(cc*cc)*(1.0/(np.exp( (hh*nu)/(kk*temp)) - 1.0))
    try:    # assume temp is an array
        bbnu = ones_like(temp)
        bbnu = 2.0 * co.h.cgs * nu**3/(co.c.cgs**2)*(1.0/(exp(x) - 1.0))
        # make sure we correct the function the same as if INT
        bbnu[x == 0] = 0.0
        ind = x < 1e-5 
        bbnu[ind] = 2.0 * nu**2 * co.k_B.cgs * temp[ind] / co.c.cgs**2
        ind = x > 20.0
        bbnu[ind] = 2.0 * co.h.cgs * nu**3 / co.c.cgs**2 * exp(-1.0*x[ind])
    except (TypeError): # if TypeError is raised, its an integer
                        # so treat it like one
        #~ print x
        if x.value == 0:
            bbnu = 0.0
        elif x.value < 1e-5:
            bbnu = 2.0 * nu**2 * co.k_B.cgs * temp / co.c**2
        elif x.value > 20.0:
            bbnu = 2.0 * co.h.cgs * nu**3 / (co.c.cgs**2) * exp(-1.0*x)
        else:
            bbnu = 2.0 * co.h.cgs * nu**3 / (co.c.cgs**2)*(1.0/(exp(x) - 1.0))
    return bbnu.to(un.Jy)

# create wavelength grid
def create_wavelength_grid(wbound = [0.01, 7., 50., 5.0e4], nwpoint = [50, 40, 40]):
    """
    
    Creates a wavlength grid with different number of points
    within differend intervals. 
    
    Input units in micro meter
    
    Output units in Hz and micro meter
    
    """
    nwav = nwpoint[0]
    # Magnus : what the hell is going on here?!?!
    # A frequency grid is created, with extra points 
    wav0  = wbound[0] * (wbound[1]/wbound[0])**(arange(nw[0], dtype=float64) / nw[0])
    for ipart in range(1,len(nw)-1): 
        dum      = wbound[ipart] * (wbound[ipart+1]/wbound[ipart])**(arange(nw[ipart], dtype=float64) / nw[ipart])
        wav = np.append(wav, dum)
    
    ipart      = len(nw)-1
    dum        = wbound[ipart] * (wbound[ipart+1]/wbound[ipart])**(arange(nw[ipart], dtype=float64) / (nw[ipart]-1.))
    wav   = np.append(wav, dum)
    # number of wavelength points
    nwav  = wav.shape[0]
    # convert wavelength in microns to frequency in Hz
    freq  = cc / (wav*1e-4)
    # number of frequency points
    nfreq = nwav
    #
    # wtf? why whould you have to sort the frequencies?
    # just flip the array, then? why this?
    # yes you could use (freq = flipud(freq) lam = flipud(wav))
    # this is about 25% faster for the whole operation
    ssubs = freq.ravel().argsort() # sort freq
    #
    freq = freq[ssubs]
    lam = wav[ssubs]
    
    return freq, lam

# read opacities

def get_opac(infile = None):
    """# Read opacities"""
    #~ opacdir = '/disks/chem9/harsono/opacities/'
    finp = os.path.join(OPACDIR, infile)
    
    try:
        data = np.genfromtxt(finp, skip_header=3, comments='#',invalid_raise=False)
    except StandardError:
        os.system('ls %s'%opacdir)
        print 'Cannot find the opacfile'
        sys.exit()

    # find the nans
    isubs = (np.isnan(data[:,0])).nonzero()[0]
    if len(isubs) > 1:
        print 'Found NaNs...'
        data = np.delete(data,isubs,0)
    
    return {'lam':data[:,0],'kabs':(data[:,1]),'kscat':(data[:,2])}

# interpolate_opacities - what it says, interpolate a given set of 
#   opactities at given wavelength to another grid of wavelengths
#   should be able to extrapolate a constant or something...

def interpolate_opacities(opac = None, freq = None):
    """
    Function grabbed from Daniel Harsono's file: transphere_helpers.py
    
    Interpolate the opacity tables to the frequency grid
    Takes one dictionary (opacities) and one array (frequencies) 
    as input.
    
    This could be optimize/simplified a bit using array operations 
    and/or built in interoplation function, but I wont to it now...
    
    """
    try:
       nfreq = freq.shape[0]
    except NameError, AttributeError:
        print ('Supply the frequencies (as an array)')
        raise(StandardError)
    
    # empty arrays to store things in...
    kabs = zeros(nfreq)
    kscat = zeros(nfreq)
    
    wl_micron = cc/freq * 1e4 # microns
    
    for ifreq in range(nfreq):        
        lampk = wl_micron[ifreq]

        # if it extends below the wavelength in the opacity file
        # set it to constant at the lowest wavenelgth in the 
        # opactiy file        
        if lampk < opac['lam'].min(): 
            kabs[ifreq] = opac['kabs'][0]
            kscat[ifreq] = opac['ksca'][0]
        # if it extends beyond the wavelength in the opacity file
        # set it to constant at the highest wavenelgth in the 
        # opactiy file        
        elif lampk > opac['lam'].max():
            kabs[ifreq] = opac['kabs'][-1]
            kscat[ifreq] = opac['ksca'][-1]
        # else it is within the interval, and we can interpolate 
        # the 
        else:
            # Find kappa
            # 
            isub = (opac['lam'] > lampk).nonzero()[0][0]-1
            eps = (log10(lampk) - log10(opac['lam'][isub])) / (log10(opac['lam'][isub+1])-log10(opac['lam'][isub]))
            # calculate the Absorption Coefficient
            kappa = (1.0-eps)*log10(opac['kabs'][isub]) + eps*log10(opac['kabs'][isub+1])
            kappa = 10.0**(kappa)
            kabs[ifreq] = kappa
            # calculate the Scattering Coefficient
            kappa = (1.0-eps)*log10(opac['ksca'][isub]) + eps*log10(opac['ksca'][isub+1])
            kappa = 10.0**(kappa)
            kscat[ifreq] = kappa
    # return the opacities as a dictionary
    # with the scattering and absorption coefficients
    return {'kabs':kabs, 'ksca':kscat}

# write opacity table
def transphere_write_opacfile(nf = -1,opac=-1):
    """
    Write the opacity input files
    """
    if nf == -1:
        print 'nf is not define here....'
        sys.exit()
        
    fout = file('dustopac.inp','w')
    fout.write('1               Format number of this file\n')
    fout.write('1               Nr of dust species\n')
    fout.write('============================================================================\n')
    fout.write('-1              Way in which this dust species is read (-1=file)\n')
    fout.write('1               Extension of name of dustopac_***.inp file\n')
    fout.write('----------------------------------------------------------------------------\n')
    fout.close()
    
    # write the dustopac_1
    
    fout = file('dustopac_1.inp','w')
    fout.write('%d  1 \n'%nf)
    fout.write('\n')
    for inu in range(nf):
        fout.write('%13.6f\n'%(opac['kabs'][inu]))
    for inu in range(nf):
        fout.write('%13.6f\n'%(0.0))
    fout.close() 

# write frequency file
def transphere_write_frequencyfile(freq):
    fout = file('frequency.inp','w')
    fout.write('%d \n'%(freq.shape[0]))
    fout.write('\n')
    for inu in range(nfreq):
        fout.write('%13.7e\n'%(freq[inu]))
    fout.close()


# find kappa at 550 nm
def find_kappa_550nm(opac);
    """
    Find kappa at 550 nm
    
    Remember : nm 1e-9, microns 1e-6 i.e. 550 nm = 0.55 microns
    
    """
    lam0 = 0.55 # in microns
    ivis = (opac['lam'] > lam0).nonzero()[0][-1]+1 # smaller than lam0
    eps = (lam0 - opac['lam'][ivis]) / (opac['lam'][ivis-1] - opac['lam'][ivis])
    kappa = (1.0*eps) * opac['kabs'][ivis-1] + eps * opac['kabs'][ivis]
    return kappa

# construct the radial grid
# this is for transphere at the moment, have to eval. different
# software, might need different points
def transphere_make_r_grid(rin=1., rout=1000., rref=100., nout=50,, nin=10):
    """
    Create logarithmically spaced grid with refinement
    If you want without refinement, just run
    logspace(log10(rin), log10(rout), num=npoints, endpoint=True)
    
    Input
    -----
    rin     : Inner radius of the grid
    rout    : Outer radius of the grid
    rref    : Transition radius between outer and inner grid
    nout    : Number of points in outer grid (outside rref)
    nin     : Number of points in inner grid (inside rref)
    
    Output
    ------
    The grid as an array of length nin + nout
    
    """
    from scipy import logspace, log10
    if any([rref<rin, rref>rout, rin>rout]):
        raise StandardError('Wrong intput, must be rin<rref<rout')
    r1 = logspace(log10(rin), log10(rref), num=nin, endpoint=False)
    r2 = logspace(log10(rref), log10(rout), num=nout, endpoint=True)
    r = concatenate( (r1, r2) )

    # Old code from Daniel, keep for reference.
    '''
    This is very complicated...
    It could be vastly simplified. Also the function needs a nref < nr
    otherwise it will create nonsens
    
    The first set of calculations just creates a logspace grid.
    Could mabye be replaced with
    r1 = logspace(log10(rin), log10(rref), num=nref, endpoint=False)
    r2 = logspace(log10(rref), log10(rout), num=nr, endpoint=True)
    r = concatenate( (r1, r2) )
 
    However these two are a bit different, for small grids 
    the differences are on the order of a few ~6%, and smaller for
    larger grids
    
    If I run
    r = transphere_make_r_grid(1,1000,50,10,100)
    
    the first 10 cells should go up to 50 AU?
    but noo...
    
    In : r[10]
    Out: 2.0733570358910702
    
    
    # Daniel's code, cannot really understand everything here...
    lgr0  = log10(rin)          # inner radius
    lgr1 = log10(rout)          # outer radius
    lgrr = log10(rref)          # refinement radius
    n1 = nref                   # refinement points
    n2 = nr-nref-1              # points left after refinement
    n = arange(nr) - n1         # 
    c = (lgr1-lgrr)/(1.0*n2)    # log spacing in the outer part part
    b = c/(lgrr-lgr0)           # log spacing in the inner part
    r = zeros(nr)               
    r[0:n1+1] = (lgrr-lgr0)*exp(b*n[0:n1+1])        # 
    r[n1+1:n1+n2+1] = (lgrr-lgr0)+c*n[n1+1:n1+n2+1]
    r = r - (r[0] + (r[0] - r[1]))
    r = r*(lgr1-lgr0)/(r.max()) + lgr0
    r = 10.0**r
    '''
    return r

# calculate the dust density given the radial grid    
def transphere_calculate_density(r, rho0, plrho, g2d=100.):
    """
    calculate the dust density as
    rho = rho0 / g2d * (r/r0)**(prlho)
        
    """
    # Correct the density profile
    rho = rho0 / g2d * (r/r0)**(prlho)
    return rho

# calculate the opacity tau
def transphere_calculate_tau(r, rho, kappa):
    """
    calculate the opacity as
    
    tau = (0.5*(rho[1:]+rho[:-1]) * kappa * (r[1:]-r[:-1])).sum()
    
    """
    
    tau = ( 0.5*(rho[1:] + rho[:-1]) * kappa * (r[1:] - r[:-1]) ).sum()
    
    return tau







