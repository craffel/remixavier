# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

'''
Runs the experiments outlined in "Estimating Timing and Channel Distortion Across Related Signals
'''

# <codecell>

import numpy as np
import estimate
import librosa
import os
import glob
import random
import scipy.io
import mir_eval
import pickle

# <codecell>

# Plot using LaTeX labels
import matplotlib as mpl
mpl.rcParams['text.usetex']=True
mpl.rcParams['text.latex.unicode']=True
mpl.rc('font',**{'size':24})

# <markdowncell>

# # Experiment 1: Artificially Corrupted TIMIT utterances
# In this experiment, a collection of 100 audio files, each consisting of 10 TIMIT utterances, is filtered and skewed randomly.  The proposed algorithm is then used to reverse the skew and filtering.

# <markdowncell>

# ## Step 1: Create the dataset
# Randomly select 10 TIMIT uttarances and concatenate them to create a single longer utterance.  Do this 100 times to create 100 examples.

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

# <markdowncell>

# ## Step 2: Apply artificial skew and filtering and estimate it
# Each utterance is randomly resampled by a factor in the range [.98, 1.02] and is filtered by a random causal filter of the form
# $$
# h[n] = \begin{cases}
# 1,& n = 0\\
# e^{-n}r[n],& 0 < n < 10\\
# 0,& n > 10
# \end{cases}
# $$
# where each $r[n] \sim \mathrm{Normal}(0, 1)$ is a Gaussian-distributed random variable with mean $0$ and variance $1$.  This resampling ratio and filter are then estimate by comparing the corrupted utterance to the original utterance.

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
    # Then set h[n] = exp[n]*Normal(0, 1) for n > 0
    random_filters[n, 1:] += np.exp( -np.arange( n_fft - 1 ) + 1 )*np.random.randn(n_fft - 1)
# Store the errors in the computed fs ratio
ratio_errors = np.zeros( len( timit_files ) )
# Store the RMS of the residual, before and after correction
residual_rms = np.zeros( len( timit_files ) )
original_rms = np.zeros( len( timit_files ) )
# Process each TIMIT file once
for n, filename in enumerate( timit_files ):
    a, fs = librosa.load( timit_files[n], sr=None )
    # First, filter it
    a_corrupted = np.convolve( a, random_filters[n] )
    # Then, skew it
    a_corrupted = librosa.resample( a_corrupted, 1, fs_ratio_true[n], 'sinc_best' )
    # Zero pad so we can compute the residual RMS now...
    a, a_corrupted = estimate.pad( a, a_corrupted )
    original_rms[n] = np.sqrt( np.mean( (a - a_corrupted)**2 ) )
    # Estimate the skew factor
    estimated_ratio = estimate.get_best_fs_ratio(a_corrupted, a, .02, 400, int(.05*a.shape[0]), int(.1*a.shape[0]))
    # Compute the relative error of this estimate
    ratio_errors[n] = np.abs( fs_ratio_true[n] - estimated_ratio)/.04
    # Reverse the skew
    a_corrupted = librosa.resample( a_corrupted, 1, 1.0/estimated_ratio )
    a, a_corrupted = estimate.pad( a, a_corrupted )
    # Compute STFTs
    A_corrupted = librosa.stft( a_corrupted, n_fft=n_fft, hop_length=n_fft/4 )
    A = librosa.stft( a, n_fft=256, hop_length=n_fft/4 )
    A = np.array(A, dtype=np.complex128)
    A_corrupted = np.array(A_corrupted, dtype=np.complex128)
    # Estimate reverse filter
    H = estimate.best_filter_coefficients( A, A_corrupted )
    # Apply reverse filter
    A_corrupted *= H
    a_corrupted = librosa.istft( A_corrupted, hop_length=n_fft/4 )
    # Compute residual and its RMS
    a, a_corrupted = estimate.pad( a, a_corrupted )
    residual = a - a_corrupted
    residual_rms[n] = np.sqrt( np.mean( residual**2 ) )
    print "Utterance {}, fs ratio error {}, relative RMS improvement {}".format(n, ratio_errors[n], original_rms[n]/residual_rms[n])

# <markdowncell>

# ## Step 3: Output results
# Plot the frequency response of the final estimated filter and print some statistics about the errors made.

# <codecell>

plt.figure( figsize=(10, 4) )
plt.plot( 1/np.abs( H[:,0] ), '.' )
plt.plot( np.abs( np.fft.rfft( random_filters[n] ) ) )
plt.ylim( [0, np.max( np.abs( np.fft.rfft( random_filters[n] ) ) ) + .1] )
plt.xlim( [0, H.shape[0]-2] )
plt.xticks( np.linspace( 0, H.shape[0], 9 ), ["${}$".format( l ) for l in np.arange( 0, 9000, 1000 )] )
plt.ylabel( '$\mathrm{Magnitude}$' )
plt.xlabel( '$\mathrm{Frequency\;(Hz)}$' )
plt.legend( ['$1/|{\\hat{H}}|$', '$|H|$'] )
plt.savefig( '../ICASSP Submission/example-filter-magnitude.pdf', bbox_inches='tight', pad_inches=.2 )

plt.figure( figsize=(10, 4) )
plt.plot( np.angle( H[:, 0] ), '.' )
plt.plot( np.angle( np.fft.rfft( random_filters[n] ) ) )
plt.ylim( [-np.pi, np.pi] )
plt.xlim( [0, H.shape[0]-2] )
plt.xticks( np.linspace( 0, H.shape[0], 9 ), ["${}$".format( l ) for l in np.arange( 0, 9000, 1000 )] )
plt.ylabel( '$\mathrm{Angle\;(Radians)}$' )
plt.xlabel( '$\mathrm{Frequency\;(Hz)}$' )
plt.legend( ['$\\angle \\hat{H}$', '$\\angle H$'] )
plt.savefig( '../ICASSP Submission/example-filter-phase.pdf', bbox_inches='tight', pad_inches=.2 )

# <codecell>

print "Average resampling error: {}%".format(np.mean(ratio_errors*100))
print "Resampling ratio recovered exaclty in {} out of {} cases".format(np.sum( ratio_errors == 0 ), ratio_errors.shape[0])
print "Average RMS without resampling/filter estimation: {}".format(np.mean(original_rms))
print "Average RMS after resampling/filter estimation: {}".format(np.mean(residual_rms))
print "Average RMS in the {} cases where resampling ratio was incorrect: {}".format(np.sum(ratio_errors != 0),
                                                                                    np.mean(residual_rms[ratio_errors != 0]))

# <markdowncell>

# # Experiment 2: Digital Music Data
# This experiment uses a collection of 10 different songs which were released on CDs alongside instrumental (no vocals) and a cappella (vocal only) versions.  The relative timing skew (possibly caused by CD recorder/reader clock drift) and channel distortion (probably caused by differences in mastering) of the instrumental and a cappella versions are estimated relative to the original mix.  The corrected versions are then subtracted from the mix and are compared against each other.

# <markdowncell>

# ## Step 1: Estimate distortion
# For each of the 10 files, the overall timing skew, local offsets, and channel distortion of the source (instrumental or a cappella mix) are estimated relative to the original mixture.  The resulting approximately channel-reversed files are written out.

# <codecell>

# Path to .wav files ripped from CDs...
data_directory = '../Dataset/Digital/'
# ... with subdirectories 1, 2, 3... each for a different song
subdirectories = [os.path.join(data_directory, d) for d in os.listdir(data_directory)]
for subdirectory in [d for d in subdirectories if os.path.isdir(d)]:
    print 'Processing file {}'.format(subdirectory)
    # Load in the full mixture
    mix, fs = librosa.load(os.path.join(subdirectory, 'M.wav'), sr=None)
    # Perform timing/channel estimation for both a cappella ('A' files) and instrumental ('I' files)
    for s in ['A', 'I']:
        # Load in the source waveform
        source, fs = librosa.load(os.path.join(subdirectory, '{}.wav'.format(s)), sr=fs)
        # Align the source to the mixture
        # Small resampling rate range because it's all digital
        mix, source_aligned = estimate.align(mix, source, fs, correlation_size=4., max_global_offset=0.,
                                             max_skew=.0001, hop=1., max_local_offset=.01)
        # Estimate the filter
        source_filtered = estimate.reverse_channel(mix, source_aligned)
        # Write out the aligned and filtered versions
        librosa.output.write_wav(os.path.join(subdirectory, '{}-aligned.wav'.format(s)), source_aligned, fs)
        librosa.output.write_wav(os.path.join(subdirectory, '{}-aligned-filtered.wav'.format(s)), source_filtered, fs)

# <markdowncell>

# ## Step 2: Compute SDR
# For each song, we compute the SDR of $m - c$, $m - c_\mathcal{R}$ and $m - c_\mathcal{F}$ relative to $s$ where we consider both the case when $c$ is the a cappella and $s$ is the instrumental and the case where $c$ is the instrumental and $s$ is the a cappella.

# <codecell>

def trim(signals, trim_length):
    '''
    First, trims each signal in signals to the minimum length.
    Then, trims of trim_length samples from the beginning and end.
    
    Input:
        signals - iterable collection of np.ndarrays
        trim_length - length in samples to trim off the beginning and end
    Output:
        signals_trimmed - tuple of np.ndarrays
    '''
    signals_trimmed = []
    min_length = min([signal.shape[0] for signal in signals])
    for signal in signals:
        signals_trimmed.append(np.array(signal[trim_length:min_length - trim_length]))
    return tuple(signals_trimmed)

# <codecell>

import glob
for subdirectory in glob.glob('../Dataset/Digital/*/'):
    for s, t in [('A', 'I'), ('I', 'A')]:
        # Cycle through c, c_R, c_F
        for variant in ['', '-aligned', '-aligned-filtered']:
            # Load in the full mixture
            mix, fs = librosa.load(os.path.join(subdirectory, 'M.wav'), sr=None)
            source, fs = librosa.load(os.path.join(subdirectory, s + variant + '.wav'), sr=None)
            mix, source = estimate.pad(mix, source)
            librosa.output.write_wav(os.path.join(subdirectory, t + variant + '-approx.wav'), mix - source, fs)

# <codecell>

# Path to .wav files ripped from CDs...
data_directory = '../Dataset/Digital/'
# ... with subdirectories 1, 2, 3... each for a different song
subdirectories = [os.path.join(data_directory, d) for d in os.listdir(data_directory)]
# Dict to store the results.  Each row has three entries, for m - c, m - c_R, m - c_F
results = {}
results['A'] = np.zeros((len(subdirectories), 3))
results['I'] = np.zeros((len(subdirectories), 3))
for row, subdirectory in enumerate([d for d in subdirectories if os.path.isdir(d)]):
    print "Processing {}".format(subdirectory)
    # Load in the full mixture
    mix, fs = librosa.load(os.path.join(subdirectory, 'M.wav'), sr=None)
    for s, t in [('A', 'I'), ('I', 'A')]:
        # Load in the target (either a cappella or instrumental)
        target, fs = librosa.load(os.path.join(subdirectory, t + '-aligned.wav'), sr=None)
        # Cycle through c, c_R, c_F
        for col, variant in enumerate(['', '-aligned', '-aligned-filtered']):
            source, fs = librosa.load(os.path.join(subdirectory, s + variant + '.wav'), sr=None)
            # Trim to same length and remove 5s from beginning and end
            mix_trimmed, source_trimmed, target_trimmed = trim([mix, source, target], 5*fs)
            # Store SDR
            results[s][row, col] = mir_eval.separation.bss_eval_sources(mix_trimmed - source_trimmed, target_trimmed)[0]
with open('digital-results.pickle', 'wb') as f:
    pickle.dump(results, f)

# <markdowncell>

# ## Step 3: Figures
# Figure 2 is an example filter response.  Figure 3 is the relative improvement of SDR due to inverting timing and channel distortion.

# <codecell>

# We'll use the second file, because it has a particularly informative filter response
subdirectory= '../Dataset/Digital/2'
# Load in mix
mix, fs = librosa.load(os.path.join( subdirectory, 'M.wav' ), sr=None)
# Use a cappella as source
source, fs = librosa.load(os.path.join( subdirectory, 'A.wav'), sr=fs)
# Make sure they are the same length
mix, source = estimate.pad( mix, source )
# Parameters
hop = fs
max_offset = fs/100
window = 4*fs
n_fft=2**13

# Estimate offset locations every "hop" samples
offset_locations, offsets, correlations = estimate.get_local_offsets( mix, source, hop, max_offset, window )

# Remove any big jumps in the offset list
offsets = estimate.remove_outliers( offsets )
# Adjust source according to these offsets
source = estimate.apply_offsets_resample( source, offset_locations, offsets )

# Make sure they are the same length again
mix, source = estimate.pad( mix, source )

# Set window size and hop size for STFTs
n_win = n_fft/2
hop = n_win/4
# Compute spectrograms
mix_spec = librosa.stft(mix, n_fft=n_fft, win_length=n_win, hop_length=hop)
source_spec = librosa.stft(source, n_fft=n_fft, win_length=n_win, hop_length=hop)
# Must be this datatype in order for minimize to work
mix_spec = np.array(mix_spec, dtype=np.complex128)
source_spec = np.array(source_spec, dtype=np.complex128)
# Compute the best filter
H = estimate.best_filter_coefficients(mix_spec, source_spec)

fig = plt.figure( figsize=(10, 4) )
host = plt.gca()
p1, = host.plot( np.abs( H ), 'k', lw=2, label='$\\mathrm{Magnitude}$' )
plt.ylim( [0, np.max( np.abs( H ) + .1) ] )
plt.xlim( [0, H.shape[0]] )
plt.xticks( np.linspace( 0, H.shape[0] - (50/44100.0)*n_fft, 12 ), ["${}$".format( l ) for l in np.arange( 0, 24, 2 )] )
plt.ylabel( '$\mathrm{Magnitude}$' )
plt.xlabel( '$\mathrm{Frequency\;(kHz)}$' )

par1 = host.twinx()
p2, = par1.plot( np.angle( H ), 'g--', lw=4, label='$\\mathrm{Phase}$' )
plt.ylim( [-np.pi, np.pi] )
plt.ylabel( '$\mathrm{Angle\;(Radians)}$' )
lines = [p1, p2]
host.legend( lines, [l.get_label() for l in lines], 'lower center' )
plt.savefig( '../ICASSP Submission/estimated-filter.pdf', bbox_inches='tight', pad_inches=.2 )

# <codecell>

with open('digital-results.pickle') as f:
    results = pickle.load(f)
# Need lots of colors for each line
colors = ['Olive', 'Sienna', 'Tomato', 'IndianRed', 'Chartreuse',
          'Maroon',  'Teal', 'b', 'SlateGray', 'SeaGreen', 'r', 'k',
          'Turquoise', 'Violet', 'Orange', 'Fuchsia', 'DarkGoldenRod',
          'Purple', 'DarkRed', 'DarkSalmon']
plt.figure( figsize=(10, 4) )
ax = plt.gca()
ax.set_color_cycle(colors)
# Plot results for acapellas, relative to column 0 which is m[n] - c[n]
plt.plot(results['A'].T - results['A'][:, 0])
# Same for instrumentals
plt.plot(results['I'].T - results['I'][:, 0])
plt.xlim( [-.01, 2.01] )
plt.xticks( [0, 1, 2], ['$m[n] - c[n]$', '$m[n] - c_\\mathcal{O}[n]$', '$m[n] - c_\\mathcal{F}[n]$'] )
plt.ylim( [-.5, 21] )
plt.yticks( [1, 7, 13, 19] )
plt.ylabel( '$\\mathrm{SDR\;Improvement}$' )
#plt.savefig( '../ICASSP Submission/SDR-improvement.pdf', bbox_inches='tight', pad_inches=.2 )
plt.legend()

# <markdowncell>

# # Experiment 3: Vinyl music data
# This experiment is similar to the one above, except that $c[n]$ is sourced from a vinyl record and both $m[n]$ and $s[n]$ come from a CD.  There are 7 examples each of vocal isolation and vocal removal.

# <markdowncell>

# ## Step 1: Estimate distortion
# For each example, we estimate the timing and channel distortion of the vinyl-sourced instrumental or a cappella relative to the CD-sourced mix.  We also perform Wiener filtering to help remove any highly nonlinear processing.

# <codecell>

data_directory = '../Dataset/Vinyl/'
# Each example has its own folder, with M.wav (mix, CD), C.wav (vinyl corrupted source), and S.wav (true source, CD)
subdirectories = [os.path.join(data_directory, d) for d in os.listdir(data_directory)]
for subdirectory in [d for d in subdirectories if os.path.isdir(d)]:
    print 'Processing file {}'.format(subdirectory)
    # Load in mixture and corrupted source
    mix, fs = librosa.load(os.path.join( subdirectory, 'M.wav' ), sr=None)
    source, fs = librosa.load(os.path.join( subdirectory, 'C.wav' ), sr=fs)
    # Align the source to the mixture
    mix, source_aligned = estimate.align(mix, source, fs, max_global_offset=0)
    # Estimate the filter
    source_filtered = estimate.reverse_channel(mix, source_aligned)
    # Write out aligned version
    librosa.output.write_wav(os.path.join(subdirectory, 'C-filtered.wav'), source_filtered, fs)
    mix, source_filtered = estimate.pad(mix, source_filtered)
    # Wiener filter the approximate separation
    enhanced = estimate.wiener_enhance( mix - source_filtered, source_aligned, 6 )
    # Write out approximation of the true source
    librosa.output.write_wav(os.path.join(subdirectory, 'S-approx.wav'), enhanced, fs)

# <markdowncell>

# ## Step 2: Compute SDRs
# As above, we compute the SDR of our approximated source against the true source.  We don't compare it against the true source with channel/timing distortion approximately removed because the vinyl distortion is much greater than for CDs.

# <codecell>

# Path to .wav files ripped from vinyl...
data_directory = '../Dataset/Vinyl/'
# Subdirectories 1-7 have C.wav an a cappella, 8-14 have instrumental
subdirectories = [os.path.join(data_directory, d) for d in os.listdir(data_directory)]
# Dict to store the results, each entry is an array of length 7
results = {}
results['removal'] = np.zeros((len(subdirectories)/2))
results['isolation'] = np.zeros((len(subdirectories)/2))
results['isolation-filtered'] = np.zeros((len(subdirectories)/2))
# Filter coefficients for vocal isolation
b, a = scipy.signal.iirfilter(2, 216./(fs/2), btype='highpass')
for subdirectory in [d for d in subdirectories if os.path.isdir(d)]:
    index = int(os.path.split(subdirectory)[1]) - 1
    print "Processing {}".format(subdirectory)
    # Load in the target (either a cappella or instrumental)
    target, fs = librosa.load(os.path.join(subdirectory, 'S.wav'), sr=None)
    approx, fs = librosa.load(os.path.join(subdirectory, 'S-approx.wav'), sr=None)
    # Trim to same length and remove 5s from beginning and end
    target_trimmed, approx_trimmed = trim([target, approx], 5*fs)
    # Figure out the task based on which file
    task = ('removal' if index < 7 else 'isolation')
    if task == 'removal':
        results[task][index] = mir_eval.separation.bss_eval_sources(approx_trimmed, target_trimmed)[0]
    elif task == 'isolation':
        results[task][index - 7] = mir_eval.separation.bss_eval_sources(approx_trimmed, target_trimmed)[0]
        approx_filtered = scipy.signal.filtfilt(b, a, approx_trimmed)
        results['isolation-filtered'][index - 7] = mir_eval.separation.bss_eval_sources(approx_filtered, target_trimmed)[0]   
with open('vinyl-results.pickle', 'wb') as f:
    pickle.dump(results, f)

# <codecell>

subdirectory = '../Dataset/Vinyl/9'
mix, fs = librosa.load(os.path.join( subdirectory, 'M.wav' ), sr=None)
source, fs = librosa.load(os.path.join( subdirectory, 'C.wav' ), sr=fs)

source = librosa.resample( source, 1, 1.010899 )
# Make sure they are the same length
mix, source = estimate.pad( mix, source )
# Parameters for local offset estimation
hop = int(fs)
max_offset = int(fs/10)
window = int(4*fs)
# Perform one iteration
n_fft=2**13

# Estimate offset locations every "hop" samples
offset_locations, offsets, correlations = estimate.get_local_offsets( mix, source, hop, max_offset, window )

# Remove any big jumps in the offset list
offsets = estimate.remove_outliers( offsets )

plt.figure( figsize=(10, 4) )
plt.imshow( np.flipud(correlations.T), aspect='auto', interpolation='nearest', cmap=plt.cm.gray )
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

with open('vinyl-results.pickle') as f:
    results = pickle.load(f)
print "Vocal Removal & {:.2f} $\pm$ {:.2f} dB\\\\".format(np.mean(results['removal']), np.std(results['removal']))
print "Vocal Isolation & {:.2f} $\pm$ {:.2f} dB\\\\".format(np.mean(results['isolation']), np.std(results['isolation']))
print "Vocal Isolation (filtered) & {:.2f} $\pm$ {:.2f} dB\\\\".format(np.mean(results['isolation-filtered']),
                                                                       np.std(results['isolation-filtered']))

