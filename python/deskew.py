# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

'''
Given two audio signals with related content, align them in time.
'''
# 2013-06-28 Dan Ellis dpwe@ee.columbia.edu + Colin Raffel craffel@gmail.com

# <codecell>

import numpy as np
import scipy.signal
import librosa

# <codecell>

def deskew(dm, dr, sr=44100):
    '''
    [y,a,b] = deskew(mix,ref,sr,PLOT)
    y is a version of mix that temporally aligns as well as possible
    with ref, which involves delaying the start by o and stretching
    the timebase by factor k.
    '''
    
    # convert to mono
    if dm.ndim > 1:
        dmm = np.sum( dm, axis=0 )
    else:
        dmm = dm
    if dr.ndim > 1:
        drm = np.sum( dr, axis=0 )
    else:
        drm = dr
    
    # This is a cut down version of skewview
    # Set the skewview parameters
    xcorrwinsec = 10.0
    xcorrhopsec = 2.0
    xcorrmaxlagsec = 2.0
    backoffsec = 0.5
    xcorrpeakthresh = 0.2
    fitthresh = 2.0
    
    # Find best global xcorr to nearest 1 ms
    n, _, _ = find_skew(dmm, drm, [], np.round(sr/1000))
    initialdelay = n/sr
    
    # apply backoff so final time offset is most likely to be a trim,
    # not insert
    n = n - np.round(backoffsec*sr)
    initialdelay = n/sr
    
    if n > 0:
        # chop off first part of mix
        dmm = dmm[n:]
    else:
        # dm starts midway through dr, pre-pad it with zeros to compensate
        dmm = np.append( np.zeros( -n ), dmm )
    
    xcorrwin = int( np.round(sr*xcorrwinsec) )
    xcorrmaxlag = int( np.round(sr*xcorrmaxlagsec) )
    xcorrhop = int( np.round(sr*xcorrhopsec) )
    
    Z, E = stxcorr(dmm, drm, xcorrwin, xcorrhop, xcorrmaxlag)
    # normalized xcorr
    ZN = Z*np.tile(1/E, (Z.shape[0], 1))
    
    zmaxpos = np.argmax( np.abs( ZN ), axis=0 )
    zmax = np.abs( ZN[zmaxpos, np.arange( ZN.shape[1] )] )
    
    # remove points where correlation is much lower than peak
    np.delete( zmaxpos, np.nonzero(zmax < xcorrpeakthresh*np.max(zmax)) )
    # actual times that corresponds to
    zmaxsec = (zmaxpos - xcorrmaxlag - 1)/(1.0*sr)
    # hence best linear fit?
    tt = np.arange( zmaxpos.shape[0] )*xcorrhop/sr
    # This doesn't ignore outliers like DAn's code does
    a, b = np.polyfit(tt, zmaxsec, 1)
    
    #scipy.io.savemat( 'workspace', {('py_' + k):v for k, v in locals().items() if isinstance(v, np.ndarray)} )
    
    a = a + 1
    
    n2 = np.round(b*sr)
    # total skew
    n = n + n2
    
    if dm.ndim > 1:
        # apply to stereo input
        if n > 0:
            dm = dm[:, :n]
        else:
            dm = np.hstack( [np.zeros( (dm.shape[0], -n) ), dm] )
        y = np.zeros( (dm.shape[0], dm.shape[1]*a) )
        for n, channel in enumerate( dm ):
            y[n, :] = scipy.signal.resample( channel, a*channel.shape[0] )
    else:
        # apply to stereo input
        if n > 0:
            dm = dm[:n]
        else:
            dm = np.append( np.zeros( -n ), dm )
        y = scipy.signal.resample( dm, a*dm.shape[0] )
    
    return y

# <codecell>

def stxcorr(X, Y=np.array([]), W=2000, H=None, N=None):
    '''
    [Z,E] = stxcorr(X,Y,W,H,N)  short-time cross-correlation
    X and Y are segmented into W-point stretches, hopped by H samples. 
    The cross-correlation of these, out to +/- N points, are returned 
    as columns of Z.
    Result is timing of X relative to Y, so that if X is equal to
    some padding samples followed by Y (i.e., Y delayed) then the
    cross correlation has peak values beyond its midpoint (in each
    column of Z).
    E returns theoretical maximum (harmonic mean energy) for each
    col of Z.
    '''

    if H is None:
        H = np.round( W/2 )
    if N is None:
        N = W

    if N > W:
        N = W
    
    npts = 2*N + 1
    LX = X.shape[0]
    if Y.shape != (0,):
        LX = np.min([LX, Y.shape[0]])
    else:
        Y = np.array( X )
    
    nfr = 1 + np.floor((LX - W)/H)
    
    Z = np.zeros((npts, nfr))
    E = np.zeros((1, nfr))
    
    nframes = 1 + int( (LX - W)/H )
    xx = np.zeros( (2*W, nframes) )
    E = np.zeros( nframes )
    for n in xrange( nframes ):
        XX = np.fft.fft( X[n*H:n*H + W], 2*W )
        YY = np.fft.fft( Y[n*H:n*H + W], 2*W )
        E[n] = np.sqrt( np.sum( np.abs( XX )**2 )*np.sum( np.abs( YY )**2 ))/(2*W)
        xx[:, n] = np.fft.ifft( XX*np.conj( YY ) ).real
    
    Z = np.vstack( [xx[-N:, :], xx[:N + 1, :]] )
    
    return Z, E

# <codecell>

def find_skew(test, ref, search_range=np.array([]), resolution=16):
    '''
    [n,xc,ll] = find_skew(test, ref, range, resolution, dosquare)
    <test> is a test waveform; <ref> is a reference we are
    comparing it to.  <n> returns the index of <test> within
    <ref>; if <test> is ref but with some samples prepended (i.e.,
    <test> is a delayed version of <ref>), <n> returns as a
    positive value, the location of the maximum of the
    cross-correlation function, or equivalently the sample index
    within <ref> that best lines up with the first sample of <test>.
    <range> (default: length(ref)/2) specifies the absolute maximum value
    to search for <n>, or, if a two-element vector, the min and max values
    to try for <n>.
    <resolution> is the required accuracy in samples.  Waveforms
    will be downsampled to do no better than this (0 = no
    downsampling, default = 16).
    <dosquare> set to 1 causes cross-correlation to be performed
    on squared signals (gives better results for signals with
    small amounts of clock drift).
    <xc> returns the actual cross-correlation function, and <ll>
    the corresponding lag indices
    '''

    if resolution > 1:
        #test = scipy.signal.resample(test, test.shape[0]/(1.0*resolution))
        test = librosa.resample( test, test.shape[0], test.shape[0]/(1.0*resolution) )
        ref = librosa.resample( ref, ref.shape[0], ref.shape[0]/(1.0*resolution) )
        search_range = np.round(search_range/resolution)
    else:
        resolution = 1;
    
    # Make signals non-negative
    test = test**2
    ref = ref**2
    
    if search_range.shape[0] == 0:
        rangemin = -ref.shape[0] + 1
        rangemax = test.shape[0]
    elif search_range.shape[0] == 1:
        rangemin = -np.abs(search_range)
        rangemax = np.abs(search_range)
    else:
        rangemin = search_range[0]
        rangemax = search_range[1]
    
    xc = np.zeros(rangemax - rangemin + 1)
    
    # cross-correlate by convolution
    xcr = scipy.signal.fftconvolve(ref[::-1], np.append( test, np.zeros(ref.shape[0]) ))[:ref.shape[0] + test.shape[0]]
    # actual lags this result corresponds to
    lagmin = -ref.shape[0] + 1
    
    xcmax = np.min(lagmin + np.nonzero( np.abs(xcr) == np.max(np.abs(xcr)))[0])
    
    # where in ll does the first value to be returned in xcvals occur?
    offset = rangemin - lagmin
    
    # which points in ll will we copy?
    touse = np.arange(np.max([0,-offset]), np.min([xc.shape[0], xcr.shape[0] - offset]))
    
    xc[touse] = xcr[offset + touse]
    ll = resolution*np.arange( rangemin, rangemax );
    
    n = resolution*xcmax
    
    return n, xc, ll

# <codecell>

if __name__=='__main__':
    a, fs = librosa.load( '../Data/wild-mix.wav', sr=None )
    b, fs = librosa.load( '../Data/wild-instr.wav', sr=None )
    c = deskew( a, b, fs )
    librosa.output.write_wav( 'wild-mix-aligned.wav', c, fs )

