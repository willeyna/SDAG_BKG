import scipy.special as spec
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import datetime
import sys
import glob
import os 
from scipy.optimize import minimize
from scipy.special import factorial
from sklearn.neighbors import KNeighborsClassifier

#resolution of healpy grid
NSIDE = 256
NPIX = hp.nside2npix(NSIDE)

#angle array of every point on the sky
m = np.arange(NPIX)
theta, phi = hp.pix2ang(nside=NSIDE, ipix=m,lonlat=True)
###################################################################

#loads in parameters for an Energy function brought in from a plot in https://arxiv.org/pdf/1306.2309.pdf
params = np.load('Params.npy')

# MISC UTILITY FUNCTIONS ---

def bisection(array,value):
    '''
    Bisecting sort algorithm for signifigance calculations
    
    Given an ``array`` , and given a ``value`` , returns an index j such that ``value`` is between array[j]
    and array[j+1]. ``array`` must be monotonic increasing. j=-1 or j=len(array) is returned
    to indicate that ``value`` is out of range below and above respectively. 
    --some dude on stack overflow, google it'''
    n = len(array)
    if (value < array[0]):
        return -1
    elif (value > array[n-1]):
        return n

    jl = 0# Initialize lower
    ju = n-1# and upper limits.
    while (ju-jl > 1):# If we are not yet done,
        jm=(ju+jl) >> 1# compute a midpoint with a bitshift
        if (value >= array[jm]):
            jl=jm# and replace either the lower limit
        else:
            ju=jm# or the upper limit, as appropriate.
        # Repeat until the test condition is satisfied.
    if (value == array[0]):# edge cases at bottom
        return 0
    elif (value == array[n-1]):# and top
        return n-1
    else:
        return jl

#used in kent dist; spherical dot product
def sph_dot(th1,th2,phi1,phi2):
    return np.sin(th1)*np.sin(th2)*np.cos(phi1-phi2) + np.cos(th1)*np.cos(th2)


def sigmoid(x, c = 3, a = .5):
    return 1 / (1 + np.exp(-c * (x - a)))

def poisson(u, x):
    y = np.exp(-u) * (u**x)/(factorial(x))
    return y

def gaussian(x,sig,mu):
    return np.exp(-np.power(x-mu,2)/(2*np.power(sig,2)))

#Takes in the track x,y,error and an array of every possible theta,phi and returns a 'signal term'.
#Uses Kent distribution/ spehrical gaussian
def evPSFd(nue,numu):
    nue=np.deg2rad(nue)
    numu=np.deg2rad(numu)

    kappa = 1./(nue[2])**2

    log_dist = np.log(kappa) - np.log(2*np.pi) - kappa + kappa*sph_dot(np.pi/2-nue[1], np.pi/2-numu[1], nue[0], numu[0])

    return np.exp(log_dist)

#takes in a number of tracks and returns an estimated error
def err(count = 1, topo = 'T'):
    if topo.upper() == 'T':
        error = .3*(np.random.rand(count)-.5)+.7
    elif topo.upper() == 'C':
        error = 6*(np.random.rand(count)-.5)+12
    else:
        print('Invalid Topology; Please use T for track error and C for cascade error')
        return
    return error

#power law function
def Efunc(E, a, b):
    return a * E**b

def pd2sig(p):
    return np.sqrt(2)*spec.erfinv(1-p)

#p-value with bisecting sort algo.
def p_value(x, bkg):
    ##bisection sorting 
    j = bisection(bkg,x)
    #all edge cases are handled inside the bisection function 
    #returns 0 if none are above as it should
    return bkg[j+1:].shape[0]/ bkg.shape[0]

def sigmoid(x, c = 3, a =.5):
    return 1 / (1 + np.exp(-c * (x - a)))

####################################################################################

# EVENT SAMPLING/ ENERGY FUNCTIONS

#Takes in a lower energy bound, power slope, and # between 0 and 1 and returns an energy
#Energy sampler created by math magic
def EventE(a,g,cdf):
    output = ((-g+1) * cdf/(((-a**(-g+1)) / (-g+1))**-1)) + a**(-g+1)
    return output**(1/(-g+1))

#fairly outdated yet now working event sampler
def ev_maker(ninj_t=0, ntrack = 200,ninj_c = 0, ncasc = 70, in_ra = 45, in_dec= 60):
    #we define 0 to be a track and 1 to be a cascade
    topos = [0, 1]

    #creates ntrack random tracks across our sky in a (ntrack, 3) array with (lon, lat, error)
    bkg_t = np.random.randint(0,high=NPIX,size=ntrack)
    bkg_t_err = err(ntrack, 'T')

    bkg_c = np.random.randint(0,high=NPIX,size=ncasc)
    #larger error for cascades
    bkg_c_err = err(ncasc, 'C')

    t_lon, t_lat = hp.pix2ang(NSIDE, bkg_t, lonlat=True)
    c_lon, c_lat = hp.pix2ang(NSIDE, bkg_c, lonlat=True)
    t_lon -= 180
    c_lon -= 180


    bkg_t = np.vstack((t_lon, t_lat, bkg_t_err)).T
    bkg_c = np.vstack((c_lon, c_lat, bkg_c_err)).T

    inj_t = np.vstack((np.random.normal(in_ra,err(),ninj_t),np.random.normal(in_dec,err(),ninj_t),err(ninj_t))).T
    inj_c = np.vstack((np.random.normal(in_ra,err(topo = 'C'),ninj_c),np.random.normal(in_dec,err(topo = 'C'),ninj_c),err(ninj_c, 'C'))).T

    #adds injected tracks to our background tracks from earlier
    tracks = np.vstack((bkg_t,inj_t))
    cascs = np.vstack((bkg_c,inj_c))
    ev = np.vstack([tracks, cascs])

    #gets energies for bkg and signal events of the correct sized array
    E = np.vstack([EventE(1e4, 3.7, np.random.random([ntrack + ncasc,1])),EventE(1e4, 2, np.random.random([ninj_t + ninj_c,1]))])

    #adds energy to events
    ev = np.concatenate([ev, E], axis = 1)

    #creates an array of 1s and 0s that shows topology for the events
    topos = np.vstack([np.zeros(shape = [ninj_t+ntrack,1]), np.ones(shape = [ninj_c+ncasc,1])])

    ev = np.concatenate([ev, topos], axis = 1)

    return ev

'''
Input: Energy in GeV, topology of neutrino (track = 1, casc = 0) , signal boolean (false = background atmospheric distribution)

Output: Percent of events that this kind of event is at the given energy
'''
#Given sig/bkg and topology, gives the percent of events of this kind at an energy
def PercentE(E, t_c, signal = True):
    perc = np.zeros_like(E)

    track_b = Efunc(E, *params[0])
    track_s = Efunc(E, *params[1])
    casc_b = Efunc(E, *params[2])
    casc_s = Efunc(E, *params[3])

    summed = (track_b + track_s + casc_b + casc_s)

    for i in range(t_c.shape[0]):
        if t_c[i]:

            if signal:
                perc[i] = casc_s[i]/ summed[i]
            else:
                perc[i] = casc_b[i]/ summed[i]
            pass

        elif not t_c[i]:
            if signal:
                perc[i] = track_s[i]/ summed[i]
            else:
                perc[i] = track_b[i]/ summed[i]
            pass

    return perc


######################################################################

# METHODS FOR SIGNAL DETECTION
# FOR EACH METHOD TO FUNCTION PROPERLY WITH THE HYPOTHESIS TESTING SCRIPTS THEY MUST RETURN TS AS RETURN VALUE [0]

#CLASSIC LLH WITHOUT ENERGY
def LLH_detector(tracks,cascades, in_ra = 45, in_dec = 60):
    evs = np.concatenate([tracks,cascades])
    nev = evs.shape[0]
    B = 1/(4*np.pi)

    S = evPSFd([evs[:,0],evs[:,1],evs[:,2]],[in_ra,in_dec])

    fun = lambda n, S, B: -np.sum(np.log( ((n/(S.shape[0]))*S) + ((1 - n/(S.shape[0]))*B) ))
    opt = minimize(fun, 10, (S,B), bounds = ((0,None),))

    n_sig = float(opt.x)
    maxllh = -float(opt.fun)
    TS = 2*(maxllh - nev*np.log(B))

    return TS, n_sig


def SMTopoAw(tracks, cascades, in_ra = 45, in_dec = 60):
    evs = np.concatenate([tracks,cascades])
    fS = PercentE(evs[:,3],evs[:,4], True)
    fB = PercentE(evs[:,3],evs[:,4], False)
    
    S = evPSFd([evs[:,0],evs[:,1],evs[:,2]],[in_ra,in_dec]) * sigmoid(fS, a = 0.5, c = 2.2)

    B = np.zeros_like(S)
    B += (1/(4*np.pi)) * fB

    fun = lambda n, S, B: -np.sum(np.log( (((n/(S.shape[0]))*S) + ((1 - n/(S.shape[0]))*B))))
    opt = minimize(fun, 10, (S,B), bounds = ((0,None),))

    injected = float(opt.x)
    maxllh = -float(opt.fun)

    TS = 2*(maxllh - np.sum(np.log(B)))

    return TS, injected

def Cascade_Prior(tracks, cascades, in_ra = 45, in_dec = 60):
    ntrack = tracks.shape[0]
    B = 1/(4*np.pi)

    S = evPSFd([tracks[:,0],tracks[:,1],tracks[:,2]],[in_ra,in_dec])

    fun = lambda n, S, B: -np.sum(np.log( ((n/(S.shape[0]))*S) + ((1 - n/(S.shape[0]))*B) ))
    opt = minimize(fun, 10, (S,B), bounds = ((0,None),))
    maxllh = -float(opt.fun)
    TS = 2*(maxllh - ntrack*np.log(B))
    #Applies a cascade prior to the traditional LLH for tracks
    PRIOR = np.sum(evPSFd([cascades[:,0],cascades[:,1],cascades[:,2]],[in_ra,in_dec]))

    TS *= PRIOR

    return TS, PRIOR

# DIFFERENT VARIENT ON CLASSIC LLH (DESCRIBED IN AN OLD POWERPOINT IN MSU ICECUBE DRIVE)
def RLLH(tracks,cascades, in_ra = 45, in_dec = 60):
    evs = np.concatenate([tracks,cascades])
    targ_ind=hp.ang2pix(NSIDE,in_ra,in_dec,lonlat=True)
    nin_ra,nin_dec=hp.pix2ang(NSIDE,targ_ind,lonlat=True)

    S = evPSFd([evs[:,0],evs[:,1],evs[:,2]],[nin_ra,nin_dec])
    B = 1/(4*np.pi)

    alpha = S > B
    ns = np.sum(alpha)
    S = S[alpha]

    TS = 2*np.sum(np.log(S/B))

    return TS, ns

# ROB'S MULTIMAP METHOD WITHOUT Energy
def MM(tracks, cascades, in_ra = 45, in_dec = 60):
    St = evPSFd([tracks[:,0],tracks[:,1],tracks[:,2]],[in_ra,in_dec])
    Sc = evPSFd([cascades[:,0],cascades[:,1],cascades[:,2]],[in_ra,in_dec])
    TS = (np.sum(St)/tracks.shape[0]) * (np.sum(Sc) / cascades.shape[0])
    return TS

# TOPOLOGY RATIO PRIOR APPLIED
# THIS VERSION DOES NOT USE knn FOR  SIGNAL COUNT; USES LLH MAXIMIZER
def TCP(tracks, cascades, in_ra = 45, in_dec = 60):
    nsc = int(round(LLH_detector(cascades)[1]))
    nst = int(round(LLH_detector(tracks)[1]))
    prior = TC[nst,nsc]

    TS0 = LLH_detector(np.concatenate([tracks,cascades]), in_ra = in_ra, in_dec = in_dec)[0]
    TS = TS0 * prior
    return TS, prior, [nst,nsc], TS0

