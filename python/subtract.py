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
    
    # Maximum offset (in samples) to search over
    offset_max = 2*fs
    # Get correlation in both directions
    mix_vs_clean = np.correlate( mix[:offset_max], source[:offset_max*2] )
    clean_vs_mix = np.correlate( source[:offset_max], mix[:offset_max*2] )
    # If the best peak was correlating mix vs clean...
    if mix_vs_clean.max() > clean_vs_mix.max():
        # Shift mix by adding zeros to the beginning of source
        offset = offset_max - np.argmax( mix_vs_clean )
        mix = np.append( np.zeros( offset ), mix )
    # Or vice versa
    else:
        offset = offset_max - np.argmax( clean_vs_mix )
        source = np.append( np.zeros( offset ), source )
    
    # Window and hop sizes
    N = 1024
    R = N/4
    
    # Make sure they are the same length again
    mix, source = pad( mix, source )
    
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
    source_filtered = librosa.istft( source_spec_filtered, n_fft=N, hop_length=R )
    # Make the same size by adding zeros
    mix, source_filtered = pad( mix, source_filtered )
    # Return remainder
    return mix - source_filtered, source_filtered
    

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
    f = 'kendrick-short'
    mix, fs = librosa.load('../Data/{}-mix-aligned.wav'.format( f ), sr=None)
    source, fs = librosa.load('../Data/{}-instr.wav'.format( f ), sr=fs)
    sep, source_filtered = separate( mix, source, fs )
    librosa.output.write_wav( '../Data/{}-sep.wav'.format( f ), sep, fs )
    enhanced = wiener_enhance( sep, source_filtered, 0 )
    librosa.output.write_wav( '../Data/{}-sep-wiener.wav'.format( f ), enhanced, fs )

