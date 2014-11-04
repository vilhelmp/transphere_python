
import errno
import os
import scipy as sp
import numpy as np
from scipy import log10, arange, array, exp, pi, ones_like, zeros
from scipy import zeros_like, size, concatenate
from ConfigParser import SafeConfigParser

# try to import AstroPy, if not raise ImportError with msg
try:
    import astropy.units as un
    import astropy.constants as co
except ImportError:
    raise ImportError('You must have the module \'astropy\' installed!')

########################################################################
# GENERAL HELP FUNCTIONS

def make_dirs(path):
    """
    # the makedirs function will raise a EEXIST error 
    # if the directory already exists, if so it returns False
    # if some other error is raise, it will raise that 
    """
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        else:
            return False
    return True

class ChangeDirectory:
    """
    Usage : 
    
    with ChangeDirectory("~/Library"):
        # we are in ~/Library
        run some code
        import subprocess
        subprocess.call("ls")
    
    after this you will still be in the previous directory
    without having to change back    
    """
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

# path to the opacities directory with the .dat file containting
# wavelength (microns), k_scattering, k_absorption
OPACDIR = "opacities"

# directory of script
__location__ = os.path.realpath(
    os.path.join(os.getcwd(), os.path.dirname(__file__)))

pc2au = 206264.806 # 1 PC in AU, or arcsecond in ???
# 206264.806 is also the number of arcsecs in a circle / 2pi

deg2rad = lambda deg : deg * pi / 180
rad2deg  = lambda rad : rad * 180 / pi


def minmax(arr):
    return [arr.min(), arr.max()]

def plancknu(nu, temp):
    """
    Calculate the planck function at a frequency and temperature
    Temperature can not be an array!
    Frequency can be an array.
    
    Input
    -----
    nu : frequency (Hz)
    temp : temperature (Kelvin) 
    
    Output
    ------
    Jansky (Jy)
    
    """
    if size(temp) > 1 :
        raise StandardError('Temperature cannot be an array')
    x =(co.h.cgs.value*nu)/(co.k_B.cgs.value*temp)
    #bbnu = 2.0*hh*nu*nu*nu/(cc*cc)*(1.0/(np.exp( (hh*nu)/(kk*temp)) - 1.0))
    try:    # assume nu is an array
        #~ bbnu = ones_like(temp)
        bbnu = 2.0 * co.h.cgs.value * nu**3 / (co.c.cgs.value**2) * (1.0/(exp(x) - 1.0))
        # make sure we correct the function the same as if INT
        bbnu[x == 0] = 0.0
        # commented out the approximations.
        # why do we need them?
        # Rayleigh-Jeans approximation (? low frequencies)
        ind = x < 1e-5
        bbnu[ind] = 2.0 * nu[ind]**2 * co.k_B.cgs.value * temp / (co.c.cgs.value**2)
        # Wiens approximation (? high frequencies)
        ind = x > 20.0
        bbnu[ind] = 2.0 * co.h.cgs.value * nu[ind]**3 / (co.c.cgs.value**2) * exp(-1.0*x[ind])
    except (TypeError): # if TypeError is raised, its an integer
                        # so treat it like one
        bbnu = 2.0 * co.h.cgs.value * nu**3 / (co.c.cgs.value**2) * (1.0/(exp(x) - 1.0))
        if x == 0:
            bbnu = 0.0
        else:
            bbnu = 2.0 * co.h.cgs.value * nu**3 / (co.c.cgs.value**2) * (1.0/(exp(x) - 1.0))
        if x == 0:
            bbnu = 0.0
        #~ # Rayleigh-Jeans approximation (? low frequencies)
        elif x < 1e-5:
            bbnu = 2.0 * nu**2 * co.k_B.cgs.value * temp / (co.c.cgs.value**2)
        # Wiens approximation (? high frequencies)
        elif x > 20.0:
            bbnu = 2.0 * co.h.cgs.value * nu**3 / (co.c.cgs.value**2) * exp(-1.0*x)
        else:
            bbnu = 2.0 * co.h.cgs.value * nu**3 / (co.c.cgs.value**2) * (1.0/(exp(x) - 1.0))
    # return Jy
    return bbnu*1e23

def plancknu_units(nu, temp):
    """
    Calculate the planck function at a frequency and temperature
    Temperature can NOT be an array.
    Frequency can be an array.
    
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
        #~ bbnu = ones_like(temp.value)
        bbnu = 2.0 * co.h.cgs * nu**3/(co.c.cgs**2)*(1.0/(exp(x) - 1.0))
        # make sure we correct the function the same as if INT
        bbnu[x == 0] = 0.0
        #~ ind = x < 1e-5 
        #~ bbnu[ind] = 2.0 * nu**2 * co.k_B.cgs * temp[ind] / co.c.cgs**2
        #~ ind = x > 20.0
        #~ bbnu[ind] = 2.0 * co.h.cgs * nu**3 / co.c.cgs**2 * exp(-1.0*x[ind])
    except (TypeError): # if TypeError is raised, its an integer
                        # so treat it like one
        if x.value == 0:
            bbnu = 0.0
        else:
            bbnu = 2.0 * co.h.cgs * nu**3 / (co.c.cgs**2)*(1.0/(exp(x) - 1.0))    
        #~ if x.value == 0:
            #~ bbnu = 0.0
        #~ elif x.value < 1e-5:
            #~ bbnu = 2.0 * nu**2 * co.k_B.cgs * temp / co.c**2
        #~ elif x.value > 20.0:
            #~ bbnu = 2.0 * co.h.cgs * nu**3 / (co.c.cgs**2) * exp(-1.0*x)
        #~ else:
            #~ bbnu = 2.0 * co.h.cgs * nu**3 / (co.c.cgs**2)*(1.0/(exp(x) - 1.0))
    return bbnu.to(un.Jy)

def rayleigh_jeans_approx(nu, temp):
    """
    Rayleigh-Jeans Approximation of the Planck Law
    Tested against plancknu above, seems legit, can reproduce 
    Wikipedia figure at least...
    """
    if size(temp) > 1 :
        raise StandardError('Temperature cannot be an array')
    x =(co.h.cgs.value*nu)/(co.k_B.cgs.value*temp)
    try:    # assume nu is an array
        # Rayleigh-Jeans approximation (? low frequencies)
        if all(x < 1e-5):
            bbnu = 2.0 * nu**2 * co.k_B.cgs.value * temp / (co.c.cgs.value**2)
        else:
            raise StandardError('Rayleigh-Jeans approximation not applicable to all elements of input')
    except (TypeError): # if TypeError is raised, its an integer
                        # so treat it like one
        # Rayleigh-Jeans approximation (? low frequencies)
        if x < 1e-5:
            bbnu = 2.0 * nu**2 * co.k_B.cgs.value * temp / (co.c.cgs.value**2)
        else:
            raise StandardError('Rayleigh-Jeans approximation not applicable to input')
    # return Jy
    return bbnu*1e23
    
def wiens_approx(nu, temp):
    """
    Wiens Approximation of the Planck Law
    Tested against plancknu above, seems legit, can reproduce 
    Wikipedia figure at least...
    """
    if size(temp) > 1 :
        raise StandardError('Temperature cannot be an array')
    x =(co.h.cgs.value*nu)/(co.k_B.cgs.value*temp)
    try:    # assume nu is an array
        # Wiens approximation (? high frequencies)
        if all(x > 20.0):
            bbnu = 2.0 * co.h.cgs.value * nu**3 / (co.c.cgs.value**2) * exp(-1.0*x)
        else:
            raise StandardError('Wiens approximation not applicable to all elements of input')
    except (TypeError): # if TypeError is raised, its an integer
                        # so treat it like one
        # Wiens approximation (? high frequencies)
        if x > 20.0:
            bbnu = 2.0 * co.h.cgs.value * nu**3 / (co.c.cgs.value**2) * exp(-1.0*x)
        else:
            raise StandardError('Wiens approximation not applicable to input')
    # return Jy
    return bbnu*1e23
    
def plancknu_tarray(nu, temp):
    """
    Calculate the Planck Law for frequency(ies) and temperature(s)
    i.e. temperature CAN be an array here
    """
    plancks = array([plancknu(nu, t) for t in temp])
    return plancks

# create wavelength grid
def create_wavelength_grid(wl_interval = [0.01, 7., 50., 5.0e4], wl_points = [50, 40, 40]):
    """
    
    Creates a wavlength grid with different number of points
    within differend intervals. 
    
    Input units in micro meter
    
    Output units in Hz and micro meter
    
    """
    #~ np = npoints[0]
    # Magnus : what the hell is going on here?!?!
    # A frequency grid is created, with extra points 
    # with some kind of spacing...
    # TODO : make this better (pythonic) and explain more...
    #~ print wl_interval
    #~ print wl_points
    wav  = wl_interval[0] * (wl_interval[1]/wl_interval[0])**(arange(wl_points[0], dtype='float64') / (wl_points[0]-1.))
    if len(wl_points)>1:
        for ipart in range(1,len(wl_points)-1): 
            dum = wl_interval[ipart] * (wl_interval[ipart+1]/wl_interval[ipart])**(arange(wl_points[ipart], dtype='float64') / (wl_points[ipart]-1.))
            wav = np.append(wav, dum)
        
        ipart = len(wl_points)-1
        dum = wl_interval[ipart] * (wl_interval[ipart+1]/wl_interval[ipart])**(arange(wl_points[ipart], dtype='float64') / (wl_points[ipart]-1.))
        wav = np.append(wav, dum)
    # number of wavelength points
    nwav  = wav.shape[0]
    # convert wavelength in microns to frequency in Hz
    # convert speed of light to micrometer / s, i.e. same length
    # unit as 'wav'
    freq  = (co.c.to(un.um/un.s)).value / wav
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
    """
    Read opacities
        
    """
    finp = os.path.join(__location__, OPACDIR, infile)
    
    try:
        data = np.genfromtxt(finp, skip_header=3, comments='#',invalid_raise=False)
    except StandardError:
        os.system('ls {0}'.format(OPACDIR))
        raise StandardError('Cannot find the opacfile')

    # find the nans
    isubs = (np.isnan(data[:,0])).nonzero()[0]
    if len(isubs) > 1:
        print 'Found NaNs...'
        data = np.delete(data,isubs,0)
    
    return {'lam':data[:,0],'kabs':(data[:,1]),'ksca':(data[:,2])}

# interpolate_opacities - what it says, interpolate a given set of 
#   opactities at given wavelength to another grid of wavelengths
#   should be able to extrapolate a constant or something...
def interpolate_opacities(opac = None, freq = None):
    """
    Function grabbed from Daniel Harsono's file: transphere_helpers.py
    
    Interpolate the opacity tables to the frequency grid
    Takes one dictionary (opacities) and one array (frequencies) 
    as input.
    
    This could probably be optimize/simplified a bit using array operations 
    and/or built in interoplation function, but I wont do it now...
    
    """
    try:
       nfreq = freq.shape[0]
    except NameError, AttributeError:
        print ('Supply the frequencies (as an array)')
        raise(StandardError)
    
    # empty arrays to store things in...
    kabs = zeros(nfreq)
    kscat = zeros(nfreq)
    
    # this will cause the wavlength to decrease with index
    # because frequency goes up with index.
    wl_micron = co.c.to(un.um/un.s).value/freq # microns
    # check...
    #~ print co.c.cgs.value / freq * 1e4 - wl_micron # microns
    # when we loop through it now, we are going from longer
    # wavelengths to shorter.
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
    
    Grabbed from Daniel Harsono's code, some modernizations 
    implemented...
    
    """
    if nf == -1:
        raise StandardError('nf is not define here....')
    
    # write the dustopac.inp file, will overwrite existing files
    fout = file('dustopac.inp','w')
    fout.write('1               Format number of this file\n')
    fout.write('1               Nr of dust species\n')
    fout.write('============================================================================\n')
    fout.write('-1              Way in which this dust species is read (-1=file)\n')
    fout.write('1               Extension of name of dustopac_***.inp file\n')
    fout.write('----------------------------------------------------------------------------\n')
    fout.close()
    
    # write the dustopac_1, will overwrite existing files
    fout = file('dustopac_1.inp','w')
    fout.write('{0}  1 \n\n'.format(int(nf))) # NB two line breaks here.
    fout.writelines(['{0:13.6f}\n'.format(i) for i in opac['kabs']])
    fout.writelines(['{0:13.6f}\n'.format(i) for i in zeros_like(opac['kabs'])])
    fout.close() 

# write frequency file
def transphere_write_frequencyfile(freq):
    fout = file('frequency.inp','w')
    fout.write('{0} \n\n'.format(int(len(freq))))
    #~ fout.write('%d \n'%(freq.shape[0]))
    #~ fout.write('\n')
    fout.writelines( ['{0:13.7f}\n'.format(f) for f in freq] )
    #~ fout.writelines( ['{0:e}\n'.format(f) for f in freq] )
    #~ for inu in range(nfreq):
        #~ fout.write('%13.7e\n'%(freq[inu]))
    fout.close()

# find kappa at X nm
def find_kappa(lam, opac_abs, lam0=550.0, localdust=False):
    """
    Find kappa at lam0 nm
    
    Remember : nm 1e-9, microns 1e-6 i.e. 550 nm = 0.55 microns
    
    Only works if 'lam' in the opac dictionary, and its increasing
    
    Just a linear interpolation between the nearest points to 0.55 nm
    could perhaps be done more robust, i.e. independent of the 
    direction of 'lam'
    """
    #~ lam0 = 0.55 # in microns
    
    if not (len(lam)==len(opac_abs)):
        raise StandardError('Input array not of same length.')
    
    if localdust:
        raise StandardError('\'localdust\' parameter not implemented yet.')
    elif not localdust:
        # first find between which 'lam' that lam0 lies
        nu_ind = (lam > lam0).nonzero()[0][0]
        # then calculate the opacity at 0.55 microns
        eps = (lam0 - lam[nu_ind-1]) / (lam[nu_ind] - lam[nu_ind-1])
        kappa = (1.-eps)*opac_abs[nu_ind-1] + eps*opac_abs[nu_ind]      
    """
    Original IDL code from Transphere release:
    ;
    ; Find the kappa at 550 nm
    ;
    if keyword_set(localdust) then begin
       openr,1,'dustopac_0.inp'
       nnr=0
       nnf=0
       idum=0
       readf,1,nnr
       readf,1,nnf,idum
       if nnr ne nr then stop
       if nnf ne nf then stop
       kap = dblarr(nf,2,nr)
       readf,1,kap
       close,1
       ivis   = 0
       fvis   = 1d4 * cc / 0.55d0 
       for inu=1,nf-1 do if((o.freq[inu]-fvis)*   $
                            (o.freq[inu-1]-fvis) lt 0.d0) then ivis=inu
       eps    = (fvis-o.freq[ivis-1]) / (o.freq[ivis]-o.freq[ivis-1])
       kappa  = transpose((1.d0-eps)*kap[ivis-1,0,*]+eps*kap[ivis,0,*],[2,0,1])
    endif

    Old Python CODE
    if localdust == 1:
        f=open('dustopac_0.inp','r')
        nr=int(f.readline().strip())
        nf,idum=readlines().strip().split()
        nf,idum=int(nf),int(idum)
        kap=np.zeros((nf,2,nr),float)
        for kk in range(0,nr):
            for jj in range(0,2):
                for ii in range(0,nf):
                    kap[ii,jj,kk]=float(f.readline().strip())
        f.close()
        ivis   = 0
        fvis   = 1e4 * nc.cc / 0.55e0 
        for inu in range(1,nf):
            if (o.freq[inu]-fvis)*(o.freq[inu-1]-fvis) < 0.e0: ivis=inu
        eps    = (fvis-o.freq[ivis-1]) / (o.freq[ivis]-o.freq[ivis-1])
        kappa  = (1.e0-eps)*kap[ivis-1,0,:]+eps*kap[ivis,0,:] # Removed transpose
      
    """
    return kappa

# construct the radial grid
# this is for transphere at the moment, have to eval. different
# software, might need different points
def make_r_grid(rin=1., rout=1000., rref=100., nout=30, nin=60):
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
    #~ r1 = logspace(log10(rin), log10(rref), num=nin, endpoint=False)
    #~ r2 = logspace(log10(rref), log10(rout), num=nout, endpoint=True)
    #~ r = concatenate( (r1, r2) )
    #~ rgrid = lambda rin, rout, nr: rin * (rout/rin)**(arange(nr)/(nr - 1.))
    print nin
    if nin>0:
        r_inner = logspace(log10(rin), log10(rref), num=nin, endpoint=False)
        r_outer = logspace(log10(rref), log10(rout), num=nout, endpoint=True)
        
        c_inner = r_inner[1:]/r_inner[:-1]
        c_outer = r_outer[1:]/r_outer[:-1]
        if c_inner[0] > c_outer[0]:
            # perhaps raise an error?
            print('You are not createing a refinement in the inner parts...')
        r = concatenate((r_inner, r_outer))
    elif not nin:
        r = logspace(log10(rin),log10(rout), num=nout, endpoint=True)
    #~ r_inner = rgrid(rin, rref, nin)
    #~ r_outer = rgrid(rref, rout, nout)
    # if no refinement 
    # THIS is a log-linear grid.
    # i.e. when ploting loglog() you get the same distance 
    # between the points.    
    #~ if nin == 0:
        #~ r = rin * (rout/rin)**(arange(nout)/(nout - 1.))
    # Old code from Daniel, keep for reference.
    #~ lgr0  = log10(rin)              # inner radius
    #~ lgr1 = log10(rout)              # outer radius
    #~ lgrr = log10(rref)              # refinement radius
    #~ n1 = nin                            # refinement points
    #~ n2 = nout - nin - 1                  # points left after refinement
    #~ n = arange(nout) - n1             # 
    #~ 
    #~ c = (log10(rout) - log10(rref))/(n2)    # log spacing in the outer part part
    #~ 
    #~ b = c / (log10(rref) - log10(rin))           # log spacing in the inner part
    #~ 
    #~ r = zeros(nout)               
    #~ 
    # the inner grid
    #~ r[0:n1+1] = (log10(rref) - log10(rin)) * exp(b*n[0:n1+1])        # 
    #~ 
    #~ # the outer grid
    #~ r[n1+1:n1+n2+1] = (log10(rref) - log10(rin)) + c * n[n1+1:n1+n2+1]
    #~ 
    #~ # change the inner cell to 
    #~ r = r - (r[0] + (r[0] - r[1]))
    #~ # ?
    #~ r = r * (log10(rout) - log10(rin))/(r.max()) + log10(rin)
    #~ # 
    #~ r = 10.0**r
    
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

def calculate_com_pos_plrho(r_grid = None, plrho = None):
	# commong exponent 
	pl = 3 - plrho
	r = ( .5*( r_grid[1:]**pl + r_grid[:-1]**pl ) )**(1 / pl)
	return r

# calculate the dust density given the radial grid    
def transphere_calculate_density(r, n0=None, r0=None, plrho=None, rho0=None, g2d=100., muh2 = 2.3):
    """
    calculate the dust density as
    rho = rho0 / g2d * (r/r0)**(plrho)
    
    Input
    -----
    r       : radial grid
    n0      : h2 number density at r0
    r0      : reference radius (same unit as 'r')
    plrho   : powerlaw exponent, including if it's negative
    g2d     : gas to dust ratio, default = 100.
    rho0    : if you want to input the dust density at r0
              you can do that
    
    
    Output
    ------
    n_h2, rho_gas, rho_dust
    
    
    """
    # Correct the density profile
    r_dependence  = (r/r0)**(plrho)
    # n0 is in nh2 / cm3
    #    g/cm3           nh2 / cm3 * amu/h2 * mass of amu(g)
    #~ self.rho0_gas = self.n0 * _cgs.MUH2 * _cgs.MP   # g * cm-3
    # all gas, MUH2 is including H2+He+Metals
    # apply the radial dependence to the gas density
    #~ self.rho_gas = self.rho0_gas * r_dependence     # g * cm-3
    # apply radial dependence to the H2 number density
    if n0 != None:
        n_h2 = n0 * r_dependence              # nh2 * cm-3
        rho_gas = n_h2 * muh2 * co.m_p.cgs.value  # g * cm-3
    elif rho0 != None:
        rho_gas = rho0 * r_dependence
        n_h2 = rho_gas / ( muh2 * co.m_p.cgs.value )
    else:
        raise StandardError('Need either \'n0\', or \'rho0\'')
    rho_dust = rho_gas / g2d  # g * cm-3
    return n_h2, rho_gas, rho_dust

# calculate the opacity tau
def transphere_calculate_tau(r, rho, kappa):
    """
    calculate the opacity as
    
    tau = (0.5*(rho[1:]+rho[:-1]) * kappa * (r[1:]-r[:-1])).sum()
    
    """
    
    tau = ( 0.5*(rho[1:] + rho[:-1]) * kappa * (r[1:] - r[:-1]) ).sum()
    
    return tau


def read_transphere_modelfile(modelfile):
    """
    Read the Transphere model file. Compatible with
    ConfigParser.SafeConfigParser.
    Input explained in the provided example file. Minimum example below
    
    Input
    -----
        [model]
        rstar =  5.97
        mstar = 1.
        tstar = 5780.  
        nin = 10
        nout = 60
        rref = 100
        rin = 35.9 
        rout = 17500.
        n0 = 4.9E8    
        plrho = -1.7
        r0 = 35.9
        opacfile = oh5_draine.dat
        wl_interval = 0.01, 7.0, 50.0, 5e4
        wl_points = 50, 40, 40

        [settings]
        plot = False
        nriter = 20
        convcrit = 1e-06
        ncst = 10
        ncex = 30
        ncnr = 1
        itypemw = 1
        idump = 1
    
    Output
    ------
    dict(), dict()
    
    Two dictionaries representing the [model] and [settings]
    sections in the input file. Some of the values in the dictionaries
    have been converted to relevant Python types.
    
    """
    # use the SafeConfigParser to read in the model file.
    parser = SafeConfigParser()
    parser.read(modelfile)
    # Store the two available sections in their own dictionaries
    model = dict(parser.items('model'))
    settings = dict(parser.items('settings'))
    # convert relevant input in model
    model['wl_interval'] = [float(i) for i in model['wl_interval'].split(',')]
    try:
        model['wl_points'] = [int(i) for i in model['wl_points'].split(',')]
    except AttributeError:
        model['wl_points'] = int(model['wl_points'])
    model['nin'] = int(model['nin'])
    model['nout'] = int(model['nout'])
    model['rref'] = float(model['rref'])
    model['rin'] = float(model['rin'])
    model['rout'] = float(model['rout'])
    if 'n0' in model.keys():
        model['n0'] = float(model['n0'])
    # if no number density was given, it is 
    # the gas mass density
    elif 'rho0' in model.keys():
        model['rh0'] = float(model['rh0'])
        # convert this to gas number density
        # using mean mass of 2.3
        muh2 = 2.3
        n0 = model['rho0'] / (muh2 * co.m_p.cgs.value)
        model['n0'] = n0
    model['r0'] = float(model['r0'])
    model['plrho'] = float(model['plrho'])
    model['tbg'] = float(model['tbg'])
    model['dpc'] = float(model['dpc'])
    if 'g2d' not in model.keys():
        print ('gas to dust ration not given, assume 100 ')
        model['g2d'] = 100.
    else:
        model['g2d'] = float(model['g2d'])
    if 'localdust' not in model.keys():
        model['localdust'] = False
    else:
        try:
            model['localdust'] = bool(eval( model['localdust'].title() ))
        except TypeError:
            raise TypeError('Parameter \'localdust\' error, evaluated to True or False only.')
    # convert relevant input in settings
    settings['plot'] = bool(settings['plot'])
    
    # return 
    return model, settings

####################
# write function
def write_opac(nf = -1,opac=-1, directory=None):
    """
    Write the opacity input files
    """
    #TODO needs update
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
    
    print 'FINISH Opacs'



# TODO need to make this better...
def read_transphereoutput(directory='test_model1'):
    from scipy import array
    filetoread = 'envstruct.dat'
    path_to_file = os.path.join(directory, filetoread)
    with open(path_to_file, 'r') as f:
        lines = f.read().split('\n')
    nr = int(lines[0])
    dat_envstruct = array([i.split() for i in lines[2:nr + 2]], dtype='float')


    filetoread = 'spectrum.dat'
    path_to_file = os.path.join(directory, filetoread)
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
    path_to_file = os.path.join(directory, filetoread)
    f = open(path_to_file, 'r')
    nn = int(f.readline().strip().split()[0])
    f.close()

    # Convergence history
    filetoread = 'convhist.dat'
    path_to_file = os.path.join(directory, filetoread)
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
    #~ if self.t_uplim:
        #~ ind = (Envstruct.temp >= self.t_uplim).nonzero()[0]
        #~ Envstruct.temp[ind] = self.t_uplim
    #~ elif self.t_flat:
        #~ ind = (self.radii <= self.t_flat).nonzero()[0]
        #ind = (self.Envstruct.temp >= self.t_uplim).nonzero()[0]
        #~ Envstruct.temp[ind] = Envstruct.temp[max(ind)]
    return Envstruct, Convhist



