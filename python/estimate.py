# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

'''
Functions for aligning and approximately removing channel distortion in signals with related content
'''

# <codecell>

import numpy as np
import scipy.optimize
import scipy.signal
import librosa

# <codecell>

def align_over_window(a, b, max_offset, correlation_size, a_center=None, b_center=None):
    '''Correlates two signals over a subset of their samples
    
    :parameters:
        - a : ndarray
            Some signal
        - b : ndarray
            Some other signal with some content in common with a
        - max_offset : int
            Maximum expected offset of the signals
        - correlation_size : int
            Number of samples to use in each correlate
        - a_center : int
            Index in a around which to center the window.  Default None, which means len(a)/2
        - b_center : int
            Index in b around which to center the window.  Default None, which means len(b)/2

    :returns:
        - offset : int
            The sample offset of a relative to b
        - a_vs_b : float
            Correlation of a against b
    '''
    # Default values for the window centers
    if a_center is None:
        a_center = a.shape[0]/2
    if b_center is None:
        b_center = b.shape[0]/2
    # Avoid array out of bounds
    if a_center - max_offset - correlation_size/2 < 0 or a_center + max_offset + correlation_size/2 > a.shape[0]:
        raise ValueError, "The window in a around a_center of size max_offset + correlation_size goes out of bounds"
    if b_center - correlation_size/2 < 0 or b_center + correlation_size/2 > b.shape[0]:
        raise ValueError, "The window in b around b_center of size max_offset + correlation_size goes out of bounds"
    # Centered on a_center, extract 2*max_offset + correlation_size samples 
    # (so offsets between -max_offset and max_offset)
    a_window = a[a_center - max_offset - correlation_size/2:a_center + max_offset + correlation_size/2]
    # From b, the sample window will only be correlation_size samples
    b_window = b[b_center - correlation_size/2:b_center + correlation_size/2]
    a_vs_b = scipy.signal.fftconvolve(a_window, b_window[::-1], 'same')[correlation_size/2:-correlation_size/2]
    # Compute offset of a relative to b
    offset = np.argmax(a_vs_b) - max_offset + (a_center - b_center)
    # Return offset and the max value of the correlation
    return offset, a_vs_b

# <codecell>

def get_best_fs_ratio(a, b, max_drift, steps, max_offset, correlation_size, center=1):
    '''
    Given two signals with components in common, tries to estimate the clock drift and offset of b vs a
    
    :parameters:
        - a : np.ndarray
            Some signal
        - b : np.ndarray
            Some other signal
        - max_drift : float
            max sample rate drift, in percent, e.g. .02 = 2% clock drift
        - steps : int
            Number of sample rates to consider, between -max_drift and max_drift
        - max_offset : int
            Maximum expected offset of the signals
        - correlation_size : int
            Number of samples to use in each correlate
        - center : float
            Ratio to deviate from - default 1
    Output:
        fs_ratio - fs ratio to make b line up well with a
    '''
    # Sample rate ratios to try
    fs_ratios = center + np.linspace(-max_drift, max_drift, steps + 1)
    # The max correlation value for each fs ratio
    corr_max = np.zeros(fs_ratios.shape)
    for n, ratio in enumerate(fs_ratios):
        # Resample b with this fs ratio
        b_resampled = librosa.resample(b, 1, ratio)
        # Compute the max correlation
        _, corr = align_over_window(a, b_resampled, max_offset, correlation_size)
        corr_max[n] = corr.max()
    # Choose ratio with the highest correlation value
    return fs_ratios[np.argmax(corr_max)]

# <codecell>

def apply_offsets_resample(b, offset_locations, offsets):
    '''
    Adjust a signal b according to local offset estimations using resampling
    
    :parameters:
        - b : np.ndarray
            Some signal
        - offset_locations : np.ndarray
            locations, in samples, of each local offset estimation
        - offsets : np.ndarray
            local offset for the corresponding sample in offset_locations
    :returns:
        - b_aligned : np.ndarray
            b with offsets applied
    ''' 
    assert offset_locations.shape[0] == offsets.shape[0]
    # Include signal boundaries in offset locations
    offset_locations = np.append(0, np.append( offset_locations, b.shape[0]-100 ))
    # Allocate output signal
    b_aligned = np.zeros(np.int(np.sum(np.diff(offset_locations)) + np.max(np.abs(offsets))))
    # Set last offset to whatever the second to last one was
    offsets = np.append(offsets, offsets[-1])
    current = 0
    # !!!!!!!!!!!!!!!!!!
    # Should zip here
    # !!!!!!!!!!!!!!!!!!
    for n, offset in enumerate(offsets):
        start = offset_locations[n]
        end = offset_locations[n + 1]
        # Compute the necessary resampling ratio to compensate for this offset
        ratio = 1 + (-offset + start - current)/(end - start)
        # Resample this portion of the signal, with some padding at the end
        resampled = librosa.resample(b[start:end + 100], 1, ratio)
        # Compute length and place the signal
        length = int(end - current - offset)
        b_aligned[current:current + length] = resampled[:length]
        current += length
    return b_aligned

# <codecell>

def apply_offsets_cola(b, offset_locations, offsets):
    '''
    Adjust a signal b according to local offset estimations using constant overlap-add
    
    :parameters:
        - b : np.ndarray
            Some signal
        - offset_locations : np.ndarray
            locations, in samples, of each local offset estimation
        - offsets : np.ndarray
            local offset for the corresponding sample in offset_locations
    :returns:
        - b_aligned : np.ndarray
            b with offsets applied
    '''
    assert offset_locations.shape[0] == offsets.shape[0]
    # Allocate output signal
    b_aligned = np.zeros(b.shape[0])
    # If first offset is not at time 0, add in original signal until first offset - should apply offset here
    if offset_locations[0] != 0:
        b_aligned[:offset_locations[0]] += np.linspace(1, 0, offset_locations[0])*b[:offset_locations[0]]
    # Include signal boundaries in offset locations
    offset_locations = np.append(0, np.append(offset_locations, b.shape[0]))
    # Add in shifted windowed signal windows
    for n in xrange(offsets.shape[0]):
        start = offset_locations[n]
        middle = offset_locations[n + 1]
        end = offset_locations[n + 2]
        if start + offsets[n] < 0:
            window = np.linspace(0, 1, middle - start)
            b_aligned[start:middle] += window*np.append(np.zeros(-(start + offsets[n])), b[:middle + offsets[n]])
        else:
            window = np.linspace(0, 1, middle - start)
            b_aligned[start:middle] += window*b[start + offsets[n]:middle + offsets[n]]
        if end + offsets[n] > b.shape[0]:
            window = np.linspace(1, 0, end - middle)
            b_aligned[middle:end] += window*np.append(b[middle + offsets[n]:], np.zeros((end + offsets[n]) - b.shape[0]))
        else:
            window = np.linspace(1, 0, end - middle)
            b_aligned[middle:end] += window*b[middle + offsets[n]:end + offsets[n]]
    # Add in final samples - should apply offset here
    window = np.linspace(0, 1, b_aligned.shape[0] - offset_locations[-2])
    b_aligned[offset_locations[-2]:] += window*b[offset_locations[-2]:b_aligned.shape[0]]
    return b_aligned

# <codecell>

def get_local_offsets(a, b, hop, max_offset, correlation_size):
    '''
    Given two signals a and b, estimate local offsets to fix timing error of b relative to a
    
    :parameters:
        - a : np.ndarray
            Some signal
        - b : np.ndarray
            Some other signal, should be the same size as a (that is, appropriately zero-padded)
        - hop : int
            Number of samples between successive offset estimations
        - max_offset : int
            Maximum expected offset of the signals
        - correlation_size : int
            Number of samples to use in each correlate
    Output:
        offset_locations - locations, in samples, of each local offset estimation
        local_offsets - Estimates the best local offset for the corresponding sample in offset_locations
    '''
    # Compute the locations where we'll estimate offsets
    offset_locations = np.arange(correlation_size + max_offset, a.shape[0] - (correlation_size + max_offset), hop)
    local_offsets = np.zeros(offset_locations.shape[0])
    correlations = np.zeros((offset_locations.shape[0], 2*max_offset ))
    for n, offset_location in enumerate(offset_locations):
        # Compute correlation
        local_offsets[n], correlations[n] = align_over_window(b, a, max_offset, correlation_size,
                                                              offset_location, offset_location)
    return offset_locations, local_offsets, correlations

# <codecell>

def best_filter_coefficients(M, R):
    '''
    Get the best vector H such that |M - HoR| is minimized, where M, H, R are complex
    
    :parameters:
        - M : np.ndarray
            STFT matrix, shape = (nbins, nframes)
        - R : np.ndarray
            STFT matrix, shape = (nbins, nframes)
    Output:
        - H : np.ndarray
            vector, shape = (nbins, 1)
    '''
    # Must be this datatype in order for minimize to work
    M = np.array(M, dtype=np.complex128)
    R = np.array(R, dtype=np.complex128)
    # 2 columns, one for real part, one for imag
    H = np.zeros((M.shape[0], 2), dtype=np.float64)
    # For each frequency bin...
    for i, (M_i, R_i) in enumerate(zip(M, R)):
        # Objective
        l1_sum = lambda H_i: np.abs(M_i.real + M_i.imag*1j - (H_i[0]*R_i.real + H_i[0]*R_i.imag*1j +
                                                              H_i[1]*R_i.real*1j - H_i[1]*R_i.imag)).sum()
        # Minimize to find filter value for this frequency bin
        H[i] = scipy.optimize.minimize(l1_sum,
                                       H[i],
                                       bounds=[(-100e100, 1e100), (-100e100, 1e100)],
                                       method='L-BFGS-B').x
    # Combine real and imaginary parts
    return (H[:, 0] + H[:, 1]*1j).reshape(-1, 1)

# <codecell>

def remove_outliers(x, median_size=13):
    '''
    Replaces any points in the vector x which lie outside of one std dev of the local median
    
    :parameters:
        - x : np.ndarray
            Input vector to clean
        - median_size : int
            Median filter size, default 13
    :returns:
        - x_cleaned : np.ndarray
            Cleaned version of x
    '''
    median_filtered = scipy.signal.medfilt(x, median_size)
    global_std = np.std(x)
    x_cleaned = x.copy()
    for n in xrange(x.shape[0]):
        if n == 0: continue
        if x_cleaned[n] < median_filtered[n] - global_std or x_cleaned[n] > median_filtered[n] + global_std:
            x_cleaned[n] = x_cleaned[n - 1]
    return x_cleaned

# <codecell>

def wiener_enhance(target, accomp, thresh=-6, transit=3, n_fft=2048):
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
    target_spec = librosa.stft(target, n_fft=n_fft, hop_length=n_fft/4)
    accomp_spec = librosa.stft(accomp, n_fft=n_fft, hop_length=n_fft/4)
    spec_ratio = librosa.logamplitude(target_spec) - librosa.logamplitude(accomp_spec)
    spec_ratio = (spec_ratio - thresh)/transit
    mask = 0.5 + 0.5*(spec_ratio/np.sqrt(1 + spec_ratio**2))
    return librosa.istft(target_spec*mask, hop_length=n_fft/4)

# <codecell>

def pad(a, b):
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
        c = np.zeros(a.shape[0])
        c[:b.shape[0]] = b
        return a, c
    if a.shape[0] < b.shape[0]:
        c = np.zeros(b.shape[0])
        c[:a.shape[0]] = a
        return c, b
    return a, b

# <codecell>

def align(a, b, fs, correlation_size=4., max_global_offset=2.,
          max_skew_offset=2., max_skew=.02, hop=.2, max_local_offset=.1):
    '''
    Aligns signal b to a by fixing global offset, finding optimal resampling rate and fixing local offsets
    
    :parameters:
        - a : np.ndarray
            Some signal
        - b : np.ndarray
            Some other signal to align to a
        - fs : int
            Sampling rate of the signals a and b
        - correlation_size : float
            Size, in seconds, of windows over which correlations will be taken, default 4
        - max_skew_offset : float
            Size, in seconds, of the maximum offset to consider for skew estimation
        - max_global_offset : float
            Length, in seconds, of the largest global offset to consider, default 2
        - max_skew : float
            Maximum percentage skew to consider, default .02
        - hop : float
            Time, in seconds, between successive offset estimations, default .2
        - max_local_offset : float
            Maximum offset in seconds for local each offset estimation, default .1

    :returns:
        - a_aligned : np.ndarray
            The signal a, unchanged except for possibly zero-padding
        - b_aligned : np.ndarray
            The signal b, aligned in time to match a
    '''
    # Make them the same length
    a, b = pad(a, b)
    # Fix any global offset
    if max_global_offset > 0:
        offset, _ = align_over_window(a, b, max_offset=int(max_global_offset*fs),
                                      correlation_size=int(correlation_size*fs))
        if offset < 0:
            a = np.append(np.zeros(-offset), a)
        elif offset > 0:
            b = np.append(np.zeros(offset), b)
    # Fix skew
    if max_skew > 0:     
        # Downsample to 2kHz for speed!  Doesn't make a big difference performance-wise
        if fs > 2000:
            fs_ds = 2000
            a_ds = librosa.resample(a, fs, fs_ds)
            b_ds =  librosa.resample(b, fs, fs_ds)
        else:
            a_ds = a.copy()
            b_ds = b.copy()
            fs_ds = fs
        if max_skew > .0001:
            # Get coarse estimate of best sampling rate
            fs_ratio = get_best_fs_ratio(a_ds, b_ds, max_skew, 200, int(max_skew_offset*fs_ds),
                                         int(correlation_size*fs_ds))
        else:
            fs_ratio = 1.
        # Get fine estimate
        fs_ratio = get_best_fs_ratio(a_ds, b_ds, .0001, 200, int(max_skew_offset*fs_ds),
                                     int(correlation_size*fs_ds), fs_ratio)
        b = librosa.resample(b, 1, fs_ratio)

    # Estimate offset locations every "hop" seconds
    offset_locations, offsets, correlations = get_local_offsets(a, b, int(fs*hop), int(fs*max_local_offset),
                                                                int(fs*correlation_size))

    # Remove any big jumps in the offset list
    offsets = remove_outliers(offsets)
    # Adjust source according to these offsets
    b = apply_offsets_resample(b, offset_locations, offsets)
    
    # Make sure they are the same length
    a, b = pad(a, b)
    
    return a, b

# <codecell>

def reverse_channel(a, b, n_fft=2**13, win_length=2**12, hop_length=2**10):
    '''
    Estimates the channel distortion in b relative to a and reverses it
    
    :parameters:
        - a : np.ndarray
            Some signal
        - b : np.ndarray
            Some other signal with channel distortion relative to a
        - n_fft : int
            Number of samples in each FFT computation, default 2**13
        - win_length : int
            Number of samples in each window, default 2**12
        - hop_length : int
            Number of samples between successive FFT computations, default 2**10
    
    :returns:
        - b_filtered : np.ndarray
            The signal b, filtered to reduce channel distortion
    '''
    # Compute spectrograms
    a_spec = librosa.stft(a, n_fft=n_fft, win_length=win_length, hop_length=hop_length)
    b_spec = librosa.stft(b, n_fft=n_fft, win_length=win_length, hop_length=hop_length)
    # Compute the best filter
    H = best_filter_coefficients(a_spec, b_spec)
    # Apply it in the frequency domain (ignoring aliasing!  Yikes)
    b_spec_filtered = H*b_spec
    # Get back to time domain
    b_filtered = librosa.istft(b_spec_filtered, win_length=win_length, hop_length=hop_length)
    return b_filtered

