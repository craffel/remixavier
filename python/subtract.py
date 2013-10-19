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
    a_vs_b = scipy.signal.fftconvolve( a, b[:b.shape[0]/2][::-1], 'valid' )
    b_vs_a = scipy.signal.fftconvolve( b, a[:a.shape[0]/2][::-1], 'valid' )
    # If the best peak was correlating mix vs clean...
    if a_vs_b.max() > b_vs_a.max():
        # Shift mix by adding zeros to the beginning of source
        b = np.append( np.zeros( np.argmax( a_vs_b ) ), b )
    # Or vice versa
    else:
        a = np.append( np.zeros( np.argmax( b_vs_a ) ), a )
    return a, b

# <codecell>

def apply_offsets_resample( b, offset_locations, offsets ):
    '''
    Adjust a signal b according to local offset estimations using resampling
    
    Input:
        b - Some signal
        offset_locations - locations, in samples, of each local offset estimation
        offsets - Estimates the best local offset for the corresponding sample in offset_locations
    Output:
        b_aligned - b with offsets applied
    ''' 
    assert offset_locations.shape[0] == offsets.shape[0]
    # Allocate output signal
    b_aligned = np.zeros( b.shape[0] + np.sum( offsets ) )
    # Include signal boundaries in offset locations
    offset_locations = np.append( 0, np.append( offset_locations, b.shape[0] ) )
    # Set last stretch amount as whatever the 
    offsets = np.append( offsets, offsets[-1] )
    current = 0
    for n, offset in enumerate( offsets ):
        start = offset_locations[n]
        end = offset_locations[n + 1]
        ratio = 1 + offset/(end - start)
        resampled = librosa.resample(b[start:end], 1, ratio)
        b_aligned[current:current + resampled.shape[0]] = resampled
        current += resampled.shape[0]
    return b_aligned

# <codecell>

def apply_offsets_cola( b, offset_locations, offsets ):
    '''
    Adjust a signal b according to local offset estimations using COLA
    
    Input:
        b - Some signal
        offset_locations - locations, in samples, of each local offset estimation
        offsets - Estimates the best local offset for the corresponding sample in offset_locations
    Output:
        b_aligned - b with offsets applied
    '''
    assert offset_locations.shape[0] == offsets.shape[0]
    # Allocate output signal
    b_aligned = np.zeros( b.shape[0] )
    # If first offset is not at time 0, add in original signal until first offset - should apply offset here
    if offset_locations[0] != 0:
        b_aligned[:offset_locations[0]] += np.linspace(1, 0, offset_locations[0])*b[:offset_locations[0]]
    # Include signal boundaries in offset locations
    offset_locations = np.append( 0, np.append( offset_locations, b.shape[0] ) )
    # Add in shifted windowed signal windows
    for n in xrange( offsets.shape[0] ):
        start = offset_locations[n]
        middle = offset_locations[n + 1]
        end = offset_locations[n + 2]
        if start + offsets[n] < 0:
            b_aligned[start:middle] += np.linspace(0, 1, middle - start)*np.append( np.zeros( -(start + offsets[n]) ), b[:middle + offsets[n]] )
        else:
            b_aligned[start:middle] += np.linspace(0, 1, middle - start)*b[start + offsets[n]:middle + offsets[n]]
        if end + offsets[n] > b.shape[0]:
            b_aligned[middle:end] += np.linspace(1, 0, end - middle)*np.append( b[middle + offsets[n]:], np.zeros( (end + offsets[n]) - b.shape[0] ) )
        else:
            b_aligned[middle:end] += np.linspace(1, 0, end - middle)*b[middle + offsets[n]:end + offsets[n]]
    # Add in final samples - should apply offset here
    b_aligned[offset_locations[-2]:] += np.linspace(0, 1, b_aligned.shape[0] - offset_locations[-2])*b[offset_locations[-2]:b_aligned.shape[0]]
    return b_aligned

# <codecell>

def get_local_offsets( a, b, hop, max_offset ):
    '''
    Given two signals a and b, estimate local offsets to fix timing error of b relative to a
    
    Input:
        a - Some signal
        b - Some other signal, should be the same size as a (that is, appropriately zero-padded)
        hop - Number of samples between successive offset estimations
        window - Maximum offset in samples for local each offset estimation
    Output:
        offset_locations - locations, in samples, of each local offset estimation
        local_offsets - Estimates the best local offset for the corresponding sample in offset_locations
    '''
    # Compute the locations where we'll estimate offsets
    offset_locations = np.arange( 2*max_offset, a.shape[0] - 2*max_offset, hop )
    local_offsets = np.zeros( offset_locations.shape[0] )
    # This will be filled in with values from a to compare against a range of b
    compare_signal = np.zeros( 4*max_offset )
    for n, i in enumerate( offset_locations ):
        # Fill in values from a - half of this is always zero
        compare_signal[:2*max_offset] = a[i - max_offset:i + max_offset]
        # Compute correlation
        correlation = scipy.signal.fftconvolve( compare_signal, b[i + 2*max_offset:i - 2*max_offset:-1], 'same' )[:2*max_offset + 1]
        # Compute this local offset
        local_offsets[n] = -(np.argmax( correlation ) - (max_offset + 1))
    return offset_locations, local_offsets

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

def remove_outliers( x ):
    '''
    Replaces any points in the vector x which lie outside of one std dev of the mean
    
    Input:
        x - Input vector to clean
    Output:
        x - Cleaned version of x
    '''
    lower_limit = np.mean( x ) - np.std( x )
    upper_limit = np.mean( x ) + np.std( x )
    for n in xrange( len(x) ):
        if n == 0: continue
        if x[n] < lower_limit or x[n] > upper_limit:
            x[n] = x[n - 1]
    return x

# <codecell>

def separate( mix, source, fs, n_iter=2, n_fft=2**13 ):
    '''
    Given a mixture signal and a source signal _which are NOT skewed_, 
    estimate the best filter and remove the source signal.
    
    Input:
        mix - signal mixture, size N
        source - source signal, size N
        fs - sampling rate
        n_iter - Number of iterations to run
        n_fft - FFT size
    Output:
        remainder - the remaining signal after removing source from mix
        source_filtered - the source, filtered by the channel estimation
    '''

    # Fix any gross timing offset
    mix, source = fix_offset( mix, source )
    # Set hop size for STFTs
    hop = n_fft/4
    # Make sure they are the same length again
    mix, source = pad( mix, source )
    for n in xrange(n_iter):
        # Fix large-scale skew error
        offset_scale = int(2*fs/(5.0*n + 1))
        offset_locations, offsets = get_local_offsets( mix, source, offset_scale, 2*offset_scale )
        # Remove any big jumps in the offset list
        offsets = remove_outliers( offsets )
        source = apply_offsets_cola( source, offset_locations, offsets )

        # Make sure they are the same length again
        mix, source = pad( mix, source )

        # Compute spectrograms
        mix_spec = librosa.stft( mix, n_fft=n_fft, hop_length=hop )
        source_spec = librosa.stft( source, n_fft=n_fft, hop_length=hop )
        # Compute the best filter
        H = best_filter_coefficients( mix_spec, source_spec )
        # Apply it in the frequency domain (ignoring aliasing!  Yikes)
        source_spec_filtered = H*source_spec
        
        # Get back to time domain
        source = librosa.istft( source_spec_filtered, n_fft=n_fft, hop_length=hop )
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
    f = 'wassup'
    mix, fs = librosa.load('../Data/{}-mix.wav'.format( f ), sr=None)
    source, fs = librosa.load('../Data/{}-instr.wav'.format( f ), sr=fs)
    sep, source_filtered = separate( mix, source, fs )
    librosa.output.write_wav( '../Data/{}-sep.wav'.format( f ), sep, fs )
    enhanced = wiener_enhance( sep, source_filtered, 0 )
    librosa.output.write_wav( '../Data/{}-sep-wiener.wav'.format( f ), enhanced, fs )

