# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

'''
Runs an experiment to test the effectiveness of the proposed source separation technique.
'''

# <codecell>

import numpy as np
import subtract
import librosa
import os
import glob
import random
import scipy.io

# <codecell>

data_directory = '../Dataset/Data/'
subdirectories = [os.path.join(data_directory, d) for d in os.listdir(data_directory)]
for subdirectory in [d for d in subdirectories if os.path.isdir(d)]:
    print subdirectory

# <codecell>

data_directory = '../Dataset/Data/'
subdirectories = [os.path.join(data_directory, d) for d in os.listdir(data_directory)]
for subdirectory in [d for d in subdirectories if os.path.isdir(d)]:
    mix, fs = librosa.load(os.path.join( subdirectory, 'M.wav' ), sr=None)
    source, fs = librosa.load(os.path.join( subdirectory, 'C.wav' ), sr=fs)
    sep, source_aligned = subtract.separate( mix, source, fs, n_iter=1 )
    librosa.output.write_wav( os.path.join( subdirectory, 'C-aligned.wav' ), source_aligned, fs )
    enhanced = subtract.wiener_enhance( sep, source_aligned, 6 )
    librosa.output.write_wav( os.path.join( subdirectory, 'S-approx.wav' ), enhanced, fs )

# <codecell>

data_directory = '../Dataset/Digital/'
subdirectories = [os.path.join(data_directory, d) for d in os.listdir(data_directory)]
for subdirectory in ['../Dataset/Digital/1', '../Dataset/Digital/2']:#[d for d in subdirectories if os.path.isdir(d)]:
    mix, fs = librosa.load(os.path.join( subdirectory, 'M.wav' ), sr=None)
    for (s, t) in [('A', 'I'), ('I', 'A')]:
        print d
        source, fs = librosa.load(os.path.join( subdirectory, '{}.wav'.format( s ) ), sr=fs)
        mix, source = subtract.pad( mix, source )
        librosa.output.write_wav( os.path.join( subdirectory, '{}-approx-naive.wav'.format( t ) ), mix - source, fs )
        mix_ds = librosa.resample( mix, fs, 8000 )
        source_ds =  librosa.resample( source, fs, 8000 )
        # Get fine estimate
        fs_ratio = subtract.get_best_fs_ratio( mix_ds, source_ds, .0001, 200, 1 )
        print fs_ratio
        source = librosa.resample( source, 1, fs_ratio )
        # Make sure they are the same length
        mix, source = subtract.pad( mix, source )
        librosa.output.write_wav( os.path.join( subdirectory, '{}-approx-skewed.wav'.format( t ) ), mix - source, fs )
        hop = fs
        max_offset = fs/100
        window = 4*fs
        # Perform one iteration
        n_fft=2**13
        mix, source_filtered, source_aligned = subtract.iteration(mix, source, hop, max_offset, window, n_fft)
        librosa.output.write_wav( os.path.join( subdirectory, '{}-approx-aligned.wav'.format( t ) ), mix - source_aligned, fs )
        librosa.output.write_wav( os.path.join( subdirectory, '{}-approx-aligned-filtered.wav'.format( t ) ), mix - source_filtered, fs )
        enhanced = subtract.wiener_enhance( mix-source_filtered, source_aligned, 0 )
        librosa.output.write_wav( os.path.join( subdirectory, '{}-approx-aligned-filtered-wiener.wav'.format( t ) ), enhanced, fs )

# <codecell>

# Get all TIMIT wav files
timit_files = glob.glob( '../Dataset/TIMIT/*.wav' )
output_dir = '../Dataset/TIMIT-concat/'
# Number of random, longer files to make
n_new_files = 100
# How many uttarances to include per file
n_utterances_per_file = 10
for n in xrange( n_new_files ):
    # Pick n_utterances_per_file TIMIT files
    random_files = random.sample( timit_files, n_utterances_per_file )
    # Concatenated utterance signal to output
    output_signal = np.array([])
    # Read in utterances and concatenate them
    for filename in random_files:
        utterance, fs = librosa.load( filename, sr=None )
        output_signal = np.append( output_signal, utterance )
    # Write librosa.output.write_wav
    librosa.output.write_wav( os.path.join( output_dir, "{}.wav".format( n ) ), output_signal, fs )

# <codecell>

# Get all TIMIT wav files
timit_files = glob.glob( '../Dataset/TIMIT-concat/*.wav' )
# Generate random resampling factors, from .98 to 1.02
fs_ratio_true = 1 + np.random.randint( -200, 201, len( timit_files ) )/10000.0
# Number of samples in each FFT window
n_fft = 256
# Generate a list of random filters to use, one for each utterance
random_filters = np.zeros( (len(timit_files), n_fft ) ) 
for n, random_filter in enumerate( random_filters ):
    # Start with an impulse...
    random_filters[n] = np.zeros( n_fft )
    random_filters[n, 0] = 1
    # Then set h[n] exp[n]*Normal(0, 1) for n > 0
    random_filters[n, 1:] += np.exp( -np.arange( n_fft - 1 ) + 1 )*np.random.randn(n_fft - 1)
# Store the errors in the computed fs ratio
ratio_errors = np.zeros( len( timit_files ) )
# Store the RMS of the residual, before and after correction
residual_rms = np.zeros( len( timit_files ) )
original_rms = np.zeros( len( timit_files ) )
# Process each TIMIT file once
for n, filename in enumerate( [timit_files[0]] ):
    a, fs = librosa.load( timit_files[n], sr=None )
    # First, filter it
    a_corrupted = np.convolve( a, random_filters[n] )
    # Then, skew it
    a_corrupted = librosa.resample( a_corrupted, 1, fs_ratio_true[n], 'sinc_best' )
    # Zero pad so we can compute the residual RMS now...
    a, a_corrupted = subtract.pad( a, a_corrupted )
    original_rms[n] = np.sqrt( np.mean( (a - a_corrupted)**2 ) )
    # Estimate the skew factor
    estimated_ratio = subtract.get_best_fs_ratio( a_corrupted, a, .02, 400 )
    # Compute the relative error of this estimate
    ratio_errors[n] = np.abs( fs_ratio_true[n] - estimated_ratio)/.04
    # Reverse the skew
    a_corrupted = librosa.resample( a_corrupted, 1, 1.0/estimated_ratio )
    a, a_corrupted = subtract.pad( a, a_corrupted )
    # Compute STFTs
    A_corrupted = librosa.stft( a_corrupted, n_fft=n_fft, hop_length=n_fft/4 )
    A = librosa.stft( a, n_fft=256, hop_length=64 )
    # Estimate reverse filter
    H = subtract.best_filter_coefficients( A, A_corrupted )
    # Apply reverse filter
    A_corrupted *= H
    a_corrupted = librosa.istft( A_corrupted, n_fft=n_fft, hop_length=n_fft/4 )
    # Compute residual and its RMS
    a, a_corrupted = subtract.pad( a, a_corrupted )
    residual = a - a_corrupted
    residual_rms[n] = np.sqrt( np.mean( residual**2 ) )
    print n, ratio_errors[n], original_rms[n]/residual_rms[n]

# <codecell>

n = 0

# <codecell>

import matplotlib as mpl
mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True
mpl.rc('font',**{'size':24})
plt.figure( figsize=(10, 4) )
plt.plot( 1/np.abs( H[:,0] ), '.' )
plt.plot( np.abs( np.fft.rfft( random_filters[n] ) ) )
plt.ylim( [0, np.max( np.abs( np.fft.rfft( random_filters[n] ) ) ) + .1] )
plt.xlim( [0, H.shape[0]-2] )
plt.xticks( np.linspace( 0, H.shape[0], 9 ), ["${}$".format( l ) for l in np.arange( 0, 9000, 1000 )] )
plt.ylabel( '$\mathrm{Magnitude}$' )
plt.xlabel( '$\mathrm{Frequency\;(Hz)}$' )
plt.legend( ['$1/|{\\hat{h}}|$', '$|h|$'] )
plt.savefig( '../ICASSP Submission/example-filter-magnitude.pdf', bbox_inches='tight', pad_inches=.2 )

plt.figure( figsize=(10, 4) )
plt.plot( np.angle( H[:, 0] ), '.' )
plt.plot( np.angle( np.fft.rfft( random_filters[n] ) ) )
plt.ylim( [-np.pi, np.pi] )
plt.xlim( [0, H.shape[0]-2] )
plt.xticks( np.linspace( 0, H.shape[0], 9 ), ["${}$".format( l ) for l in np.arange( 0, 9000, 1000 )] )
plt.ylabel( '$\mathrm{Angle\;(Radians)}$' )
plt.xlabel( '$\mathrm{Frequency\;(Hz)}$' )
plt.legend( ['$\\angle \\hat{h}$', '$\\angle h$'] )
plt.savefig( '../ICASSP Submission/example-filter-phase.pdf', bbox_inches='tight', pad_inches=.2 )

# <codecell>

np.save( 'ratio_errors.npy', ratio_errors )
np.save( 'original_rms.npy', original_rms )
np.save( 'residual_rms.npy', residual_rms )

# <codecell>

print np.mean( ratio_errors ), np.std( ratio_errors )
print np.sum( ratio_errors != 0 ), ratio_errors.shape[0]
print np.mean( ratio_errors[ratio_errors != 0] ), np.std( ratio_errors[ratio_errors != 0] )
print np.mean( original_rms ), np.std( original_rms )
print np.mean( residual_rms ), np.std( residual_rms )
print np.mean( original_rms[ratio_errors == 0] ), np.std( original_rms[ratio_errors == 0] )
print np.mean( residual_rms[ratio_errors == 0] ), np.std( residual_rms[ratio_errors == 0] )
print np.mean( original_rms[ratio_errors != 0] ), np.std( original_rms[ratio_errors != 0] )
print np.mean( residual_rms[ratio_errors != 0] ), np.std( residual_rms[ratio_errors != 0] )

# <codecell>

digital_results = scipy.io.loadmat( 'bss_eval/digital-results.mat' )['results']
colors = ['Olive', 'Sienna', 'Tomato', 'IndianRed', 'Chartreuse','Maroon',  'Teal', 'b', 'SlateGray', 'SeaGreen', 'r', 'k',  'Turquoise', 'Violet', 'Orange', 'Fuchsia', 'DarkGoldenRod', 'Purple', 'DarkRed', 'DarkSalmon']
plt.figure( figsize=(10, 4) )
ax = plt.gca()
ax.set_color_cycle(colors)
plt.plot( digital_results[0, :-1, :] - digital_results[0, 0, :] )
plt.plot( digital_results[1, :-1, :] - digital_results[1, 0, :] )
plt.xlim( [-.01, 2.01] )
plt.xticks( [0, 1, 2], ['$m[n] - c[n]$', '$m[n] - c_\\mathcal{R}[n]$', '$m[n] - c_\\mathcal{F}[n]$'] )
plt.ylim( [-.2, 8] )
plt.yticks( [1, 3, 5, 7] )
plt.ylabel( '$\\mathrm{SDR\;Improvement}$' )
plt.savefig( '../ICASSP Submission/SDR-improvement.pdf', bbox_inches='tight', pad_inches=.2 )

# <codecell>

subdirectory= '../Dataset/Digital/2'
mix, fs = librosa.load(os.path.join( subdirectory, 'M.wav' ), sr=None)
s = 'A'
t = 'I'
source, fs = librosa.load(os.path.join( subdirectory, '{}.wav'.format( s ) ), sr=fs)
# Make sure they are the same length
mix, source = subtract.pad( mix, source )
hop = fs
max_offset = fs/100
window = 4*fs
# Perform one iteration
n_fft=2**13

# Estimate offset locations every "hop" samples
offset_locations, offsets, correlations = subtract.get_local_offsets( mix, source, hop, max_offset, window )

# Remove any big jumps in the offset list
offsets = subtract.remove_outliers( offsets )
# Adjust source according to these offsets
source = subtract.apply_offsets_resample( source, offset_locations, offsets )

# Make sure they are the same length again
mix, source = subtract.pad( mix, source )

# Set window size and hop size for STFTs
n_win = n_fft/2
hop = n_win/4
# Compute spectrograms
mix_spec = librosa.stft( mix, n_fft=n_fft, hann_w=n_win, hop_length=hop )
source_spec = librosa.stft( source, n_fft=n_fft, hann_w=n_win, hop_length=hop )
# Compute the best filter
H = subtract.best_filter_coefficients( mix_spec, source_spec )

# <codecell>

fig = plt.figure( figsize=(10, 4) )
host = plt.gca()
p1, = host.plot( np.abs( H ), 'k', lw=2, label='$\\mathrm{Magnitude}$' )
plt.ylim( [0, np.max( np.abs( H ) + .1) ] )
plt.xlim( [0, H.shape[0]] )
plt.xticks( np.linspace( 0, H.shape[0] - (50/44100.0)*n_fft, 12 ), ["${}$".format( l ) for l in np.arange( 0, 24, 2 )] )
plt.ylabel( '$\mathrm{Magnitude}$' )
plt.xlabel( '$\mathrm{Frequency\;(kHz)}$' )

#plt.figure( figsize=(10, 4) )
par1 = host.twinx()
p2, = par1.plot( np.angle( H ), 'g--', lw=4, label='$\\mathrm{Phase}$' )
plt.ylim( [-np.pi, np.pi] )
plt.ylabel( '$\mathrm{Angle\;(Radians)}$' )
lines = [p1, p2]
host.legend( lines, [l.get_label() for l in lines], 'lower center' )
plt.savefig( '../ICASSP Submission/estimated-filter.pdf', bbox_inches='tight', pad_inches=.2 )

# <codecell>

subdirectory = '../Dataset/Data/9'
mix, fs = librosa.load(os.path.join( subdirectory, 'M.wav' ), sr=None)
source, fs = librosa.load(os.path.join( subdirectory, 'C.wav' ), sr=fs)

source = librosa.resample( source, 1, 1.010899 )
# Make sure they are the same length
mix, source = subtract.pad( mix, source )
# Parameters for local offset estimation
hop = int(fs/(5.0*n + 1))
max_offset = fs/10
window = int(4*fs/(2.0*n + 1))
# Perform one iteration
n_fft=2**13

# Estimate offset locations every "hop" samples
offset_locations, offsets, correlations = subtract.get_local_offsets( mix, source, hop, max_offset, window )

# Remove any big jumps in the offset list
offsets = subtract.remove_outliers( offsets )

plt.figure( figsize=(10, 4) )
plt.imshow( correlations.T, aspect='auto', interpolation='nearest', cmap=plt.cm.gray )
plt.plot( -offsets + max_offset + 1, 'o', mfc='none' )
offsets_to_plot = -offsets + max_offset - 100
center = np.mean( offsets_to_plot )
spread = np.max( offsets_to_plot ) - np.min( offsets_to_plot )
plt.xlim(-.5, offsets.shape[0]-.5 )
plt.ylim( center - .7*spread, center + .7*spread )
plt.xticks( np.linspace( 0, 180, 10 ) )
plt.xlabel( '$\\mathrm{Time\;(seconds)}$' )
plt.yticks( np.linspace( 2300, 5700, 6 ), ["${}$".format( int(d) ) for d in 1000*(np.linspace( 2300, 5700, 6 ) - 4410.0)/fs] )
plt.ylabel( '$\\mathrm{Offset\;(milliseconds)}$' )
plt.savefig( '../ICASSP Submission/correlation.pdf', bbox_inches='tight', pad_inches=.2 )

# <codecell>

SDRs_vinyl_removal = scipy.io.loadmat( 'bss_eval/vinyl-removal.mat' )['SDRs_vinyl_removal']
SDRs_vinyl_isolation = scipy.io.loadmat( 'bss_eval/vinyl-isolation.mat' )['SDRs_vinyl_isolation']
SDRs_vinyl_isolation_filtered = scipy.io.loadmat( 'bss_eval/vinyl-isolation-filtered.mat' )['SDRs_vinyl_isolation_filtered']
print "Removal & {:.2f} $\pm$ {:.2f} dB\\\\".format( np.mean( SDRs_vinyl_removal ), np.std( SDRs_vinyl_removal ) )
print "Isolation & {:.2f} $\pm$ {:.2f} dB\\\\".format( np.mean( SDRs_vinyl_isolation ), np.std( SDRs_vinyl_isolation ) )
print "Isolation (filtered) & {:.2f} $\pm$ {:.2f} dB\\\\".format( np.mean( SDRs_vinyl_isolation_filtered ), np.std( SDRs_vinyl_isolation_filtered ) )

