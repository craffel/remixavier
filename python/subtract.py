# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

'''
Functions for removing one signal from the other in a smart way
'''
# 2013-06-28 Dan Ellis dpwe@ee.columbia.edu + Colin Raffel craffel@gmail.com

# <codecell>

import numpy as np
import scipy.optimize
import scipy.signal
import librosa

# <codecell>

def fix_offset( a, b ):
    '''
    Given two signals with components in common, produces signals such that the offset is approximately zero

    Input:
        a - Some signal
        b - Some other signal
    Output:
        a_offset - Version of "a" with alignment fixed
        b_offset - Version of "b" with alignment fixed
    '''
    # Get correlation in both directions
    a_vs_b = scipy.signal.fftconvolve( a[::2], b[:b.shape[0]/2][::-2], 'valid' )
    b_vs_a = scipy.signal.fftconvolve( b[::2], a[:a.shape[0]/2][::-2], 'valid' )
    # If the best peak was correlating mix vs clean...
    if a_vs_b.max() > b_vs_a.max():
        # Shift mix by adding zeros to the beginning of source
        b = np.append( np.zeros( np.argmax( a_vs_b )*2 ), b )
    # Or vice versa
    else:
        a = np.append( np.zeros( np.argmax( b_vs_a )*2 ), a )
    return a, b

# <codecell>

def fix_skew( a, b, hop, max_offset ):
    '''
    Given two signals a and b, estimate local offsets to fix timing error of b relative to a
    
    Input:
        a - Some signal
        b - Some other signal, should be the same size as a (that is, appropriately zero-padded)
        hop - Number of samples between successive offset estimations
        window - Maximum offset in samples for local each offset estimation
    Output:
        a - Signal a, unchanged
        b_aligned - Signal b, with time corrections to align to a
    '''
    # Store the magnitude of local differences in this
    diff_mags = np.zeros( 2*max_offset)
    # Estimate first offset at the first hop location
    current_offset = 2*max_offset + hop
    # Allocate aligned version of b (to be filled in as we go)
    b_aligned = np.zeros( b.shape[0] )
    # Fill in the beginning
    b_aligned[:current_offset] = b[:current_offset]
    b_aligned[current_offset - hop:current_offset] *= np.arange( hop, 0, -1 )
    # Window for fading in/out the corrected portions
    fade_window = np.append(np.arange( hop ), np.arange( hop, 0, -1 ))
    # Estimate offsets until we run out of samples
    while current_offset + 2*max_offset + hop < b.shape[0]:
        # Grab current frame in the reference window
        current_range = r_[current_offset - max_offset:current_offset + max_offset]
        current_frame = a[current_range]
        for n, offset in enumerate( xrange( -max_offset, max_offset ) ):
            shifted_range = r_[current_offset + offset - max_offset:current_offset + offset + max_offset]
            diff_mags[n] = np.sum( np.abs( current_frame - b[shifted_range] ) )
        # Correct the offset
        offset = np.argmin( diff_mags ) - max_offset
        shifted_range = r_[current_offset + offset - hop:current_offset + offset + hop]
        b_aligned[current_offset - hop:current_offset + hop] += fade_window*b[shifted_range]
        current_offset += hop
        #print current_offset/float( b.shape[0] ), offset
    return a, b_aligned

# <codecell>

def best_filter_coefficients( M, R ):
    '''
    Get the best vector H such that |M - HoR| is minimized, where M, H, R are complex
    
    Input:
        M - matrix, size nbins x nframes
        R - matrix, size nbins x nframes
    Output:
        H - vector, size nbins x 1
    '''
    H = np.zeros( (M.shape[0], 2) )
    # For each frequency bin...
    for i, (M_i, R_i) in enumerate( zip( M, R ) ):
        l1_sum = lambda H_i: np.sum( np.abs( M_i.real + M_i.imag*1j - (H_i[0]*R_i.real + H_i[0]*R_i.imag*1j + H_i[1]*R_i.real*1j - H_i[1]*R_i.imag) ) )
        H[i] = scipy.optimize.minimize( l1_sum, H[i], bounds=[(-100e100, 1e100), (-100e100, 1e100)], method='L-BFGS-B' ).x
    return (H[:, 0] + H[:, 1]*1j).reshape( -1, 1 )

# <codecell>

def separate( mix, source, fs ):
    '''
    Given a mixture signal and a source signal _which are NOT skewed_, 
    estimate the best filter and remove the source signal.
    
    Input:
        mix - signal mixture, size N
        source - source signal, size N
        fs - sampling rate
    Output:
        remainder - the remaining signal after removing source from mix
        source_filtered - the source, filtered by the channel estimation
    '''
    
    # Fix any gross timing offset
    mix, source = fix_offset( mix, source )
    # Make sure they are the same length again
    mix, source = pad( mix, source )
    
    # Window and hop sizes
    N = 1024
    R = N/4
    
    # Compute spectrograms
    mix_spec = librosa.stft( mix, n_fft=N, hop_length=R )
    source_spec = librosa.stft( source, n_fft=N, hop_length=R )
    
    # Compute the best filter
    H = best_filter_coefficients( mix_spec, source_spec )
    
    # Apply it in the frequency domain (ignoring aliasing!  Yikes)
    source_spec_filtered = H*source_spec
    
    '''# Compute a highpass filter to get rid of the low end
    low_h = scipy.signal.firwin(N + 1, 200.0/fs, pass_zero=False)[:-1]
    low_H = np.fft.rfft( low_h ).reshape( -1, 1 )
    # Apply filter
    mix_spec = low_H*mix_spec
    source_spec_filtered = low_H*source_spec'''
    
    # Get back to time domain
    source = librosa.istft( source_spec_filtered, n_fft=N, hop_length=R )
    # Make the same size by adding zeros
    mix, source = pad( mix, source )
    
    # Now, fix the skew (parameters set arbitrarily)
    mix, source = fix_skew( mix, source, int(fs*.1), int(fs*.1) )
    
    # Compute spectrograms
    mix_spec = librosa.stft( mix, n_fft=N, hop_length=R )
    source_spec = librosa.stft( source, n_fft=N, hop_length=R )

    # Compute the best filter
    H = best_filter_coefficients( mix_spec, source_spec )
    
    # Apply it in the frequency domain (ignoring aliasing!  Yikes)
    source_spec_filtered = H*source_spec

    # Get back to time domain
    source = librosa.istft( source_spec_filtered, n_fft=N, hop_length=R )
    # Make the same size by adding zeros
    mix, source = pad( mix, source )
    
    # Return remainder
    return mix - source, source
    

# <codecell>

def wiener_enhance( target, accomp, thresh=-6, transit=3, n_fft=2048 ):
    '''
    Given a noisy signal and a signal which approximates the noise, try to remove the noise.
    
    Input:
        target - Noisy signal
        accomp - Approximate noise
        thresh - Sigmoid threshold, default -6
        tranist - Sigmoid transition, default 3
        n_fft - FFT length, default 2048 (hop is always n_fft/4)
    Output:
        filtered - Target, Wiener filtered to try to remove noise
    '''
    target_spec = librosa.stft( target, n_fft=n_fft, hop_length=n_fft/4 )
    accomp_spec = librosa.stft( accomp, n_fft=n_fft, hop_length=n_fft/4 )
    spec_ratio = librosa.logamplitude( target_spec ) - librosa.logamplitude( accomp_spec )
    spec_ratio = (spec_ratio - thresh)/transit
    mask = 0.5 + 0.5*(spec_ratio/np.sqrt(1 + spec_ratio**2))
    return librosa.istft( target_spec*mask, n_fft=n_fft, hop_length=n_fft/4 )

# <codecell>

def pad( a, b ):
    '''
    Given two vectors, pad the shorter one with zeros (at the end) so that they are the same size
    
    Input:
        a - vector
        b - vector
    Output:
        a - vector, padded if needed
        b - vector, padded if needed
    '''
    if a.shape[0] > b.shape[0]:
        c = np.zeros( a.shape[0] )
        c[:b.shape[0]] = b
        return a, c
    if a.shape[0] < b.shape[0]:
        c = np.zeros( b.shape[0] )
        c[:a.shape[0]] = a
        return c, b
    return a, b

# <codecell>

if __name__ == '__main__':
    # 2013-06-28 Dan Ellis dpwe@ee.columbia,edu + Colin Raffel craffel@gmail.com
    f = 'mc-paul'
    mix, fs = librosa.load('../Data/{}-mix.wav'.format( f ), sr=None)
    source, fs = librosa.load('../Data/{}-instr.wav'.format( f ), sr=fs)
    mix = mix[20*fs:40*fs]
    source = source[20*fs:40*fs]
    sep, source_filtered = separate( mix, source, fs )
    librosa.output.write_wav( '../Data/{}-sep.wav'.format( f ), sep, fs )
    enhanced = wiener_enhance( sep, source_filtered, 0 )
    librosa.output.write_wav( '../Data/{}-sep-wiener.wav'.format( f ), enhanced, fs )

