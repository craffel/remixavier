%% REMIXAVIER - Tools for recombining different mixes of a track
%
% 2013-06-28 Dan Ellis dpwe@ee.columbia.edu
%
% The Remixavier ("remix savior") project is concerned with recovering 
% the "difference" between different mixes of the same track.  For
% instance, given a full mix and an instrumental, we can try to
% recover the vocals, or given the full mix and an a cappella
% version, we can try to produce an instrumental version.  In the
% process, we can identify the precise temporal alignment between
% the two versions, which may be useful in its own right.
%
% Assuming we have a full mix M and an instrumental version I,
% under ideal conditions we could recover the vocal line V as M -
% I.  However, there are very often timing offsets and small sampling
% rate differences (clock drift) that will defeat the simple
% approach.  We estimate these timing differences with short-time 
% cross correlation (in <deskew.m>), and trim and resample to
% correct it to within a few parts per million (milliseconds of
% drift over the duration of a typical track).
% 
% But even with a perfect or near-perfect time alignment, there may
% be differences in gain, or more generally the channel frequency
% response, that will still make simple subtraction inadequate.
% Instead, we estimate the optimal equalization filter H to
% minimize the energy of M - H.I.  This is done within 
% <find_in_mix.m>, which calls <decomp_lin_win.m> to break the 
% pair of signals into short chunks (e.g. 8 second chunks every 4
% second), estimate the best coupling impulse response of each
% chunk in <decomp_lin.m>, then overlap-add the canceled residuals 
% to produce the desired difference.  This actually works by
% estimating a whitening filter for I so that the cross-correlation of
% the whitened versions of M and I is simply the coupling impulse
% response.  
%
% This approach to cancelation is inspired by the 
% <http://bass-db.gforge.inria.fr/bss_eval/ BSS_EVAL>
% procedure of Fevotte, Gribonval, and Vincent.  Essentially, we
% are finding the difference as the "artefacts residual" when M is
% considered an imperfect estimate of I.  

%% Example 1: Significant time skew and channel difference
%
% This example consists of an original instrumental track,
% digitized from a vinyl LP release, and a rap that uses the track
% as backing, taken directly from a CD.  Thus, the different signal
% paths mean that the timing is significantly different (clock
% drift of 0.1%), and the overall spectrum is very different too. 

% Load in mix and acapella as mono files
% These tracks diverge at the end (different edits), so just work
% on the first minute 
sr = 44100;
[dmix,sr] = mp3read('../Data/mc-paul-mix.mp3',[0 60*sr],1);
[dins,sr] = mp3read('../Data/mc-paul-instr.mp3',[0 60*sr],1);

% Attempt to trim and resample the full version to line up as well
% as possible with the acapella
doplot = 1;
dmr = deskew(dmix, dins, sr, doplot);
axis([0.5 55.5 1 1.5])
% It gets better when you repeat it
dmr = deskew(dmr, dins, sr);
% resampling can't handle ratios below 30 ppm, will just skip
% beyond that.

%% Channel estimation for example 1

% Do the short-time coupling filter estimation
tic; [resid, targ, filt, SNR, del, filts] = ...
      find_in_mix(dmr,dins,sr,0.013,0.003); toc

% Listen to the residual (vocals)
% (play the second 20 seconds)
ix = 20*sr+[1:20*sr];
soundsc(resid(ix,:), sr);

% Plot the time-local coupling filters (right channel)
% filter IR time base
tt = [1:size(filts,1)]/sr;
% times of individual short-time window (every 4 sec)
tw = 4.0*[1:size(filts,2)];
% plot
imagesc(tt,tw,filts'); axis('xy');
xlabel('time / sec')
ylabel('window time / sec')
title('local coupling filter impulse responses (cap -> mix)')
% scale down impulse response extremes
caxis([-2 2])

%% Example 2: Recovering instrumental, in stereo
%
% The Duffy track has the vocals in stereo, we can cancel left and
% right separately to good effect

[dmix,sr] = mp3read('../Data/Duffy.WarwickAvenue.mp3');
[dcap,sr] = mp3read('../Data/duffy_-_warwick_avenue_acapella.mp3');

% Deskew will process stereo files.  Skew is estimated from an
% internally-generated mono mix
dmr = deskew(dmix, dcap, sr);
dmr = deskew(dmr, dcap, sr);

clear resid targ

for i = 1:size(dmr,2)
  tic; [resid(:,i), targ(:,i), filt, SNR, del, filts] = ...
                find_in_mix(dmr(:,i),dcap(:,i),sr,0.006,0.003); toc
end
  
soundsc(resid(ix,:), sr);


%% Example 3: Perfectly-aligned signals, and Weiner enhancement

% Message In A Bottle is an ideal case - plain subtraction of mix
% and instrumental yeilds clean vocals.  But how does estimation do?

% Load tracks as mono
sr = 44100;
dmix = mean(wavread('../Data/message-in-a-bottle-mix.wav'),2);
dins = mean(wavread('../Data/message-in-a-bottle-ins.wav'),2);

% They shoyld be perfectly aligned already, but run deskew just in case
dmr = deskew(dmix, dins, sr);

tic; [resid, targ, filt, SNR, del, filts] = ...
      find_in_mix(dmr,dins,sr,0.013,0.003); toc

soundsc(resid(ix,:), sr);

% We can apply a "wiener filter" (scaling of spectrogram magnitude
% cells) to further reduce residual artifacts.  In particular, we
% can suppress cells where the energy in the estimated vocals is
% significantly lower than the energy in the instrumental line
% projected into the mix.  weinerenhace takes a threshold so that
% energy in the residual that is below -6 dB when compared to the 
% accompaniment is suppressed

reswf = weinerenhance(resid, targ, -6.0);
soundsc(reswf(ix,:), sr);

% We can measure SNR by canceling against the true vocals, which
% are simply the difference of dmix and dins (for this perfect example)
dvox = dmix - dins;
soundsc(dvox(ix,:), sr);  % Yes, sounds clean
[r2, t2, f2, S2, d2, fs2] = find_in_mix(resid,dvox,sr,0.010,0.003); 
%Delay = 0.000000 s
%SNR = 19.9197 dB  <-- this is our estimate of SDR
[r2, t2, f2, S2, d2, fs2] = find_in_mix(reswf,dvox,sr,0.010,0.003);
%Delay = 0.000000 s
%SNR = 16.1096 dB
% Weiner filtering introduces more artifact energy than it removes
% interference.

%% Example 4: Imogen Heap Instrumental Version
%
% This album was released with two versions of every track - a full
% mix, and an instrumental version.  Since they are derived from
% the same digital masters, there is no clock drift, although they
% are not perfectly aligned in time.  However, because each short
% segment reflects the same timing alignment, we can average the
% estimated coupling filters to further stabilize the estimation.
% Because there are a few outlier frames (degenerate estimates from
% when the vocal track is near silent), we combine across filters
% with a median instead of a mean.

[dmix,sr] = audioread('../Data/10-Aha_.m4a');
[dins,sr] = audioread('../Data/23-Aha_Instrumental_Version_.m4a');

% Deskew once just to remove any gross timing offset
dmr = deskew(dmix, dins, sr);

clear resid targ filts

% Align each channel, and store all the individual filters
for i = 1:size(dmr,2)
  tic; [resid(:,i), targ(:,i), filt, SNR, del, filts{i}] = ...
                find_in_mix(dmr(:,i),dins(:,i),sr,0.006,0.003); toc
end

soundsc(resid(ix,:), sr);

% but form a grand average filter for each side
f1 = median(filts{1}');
f2 = median(filts{2}');

% The estimated filter has a pre-echo, so trim that from the convolution
[vv,xx] = max(abs(f1));
% xx is the index of the peak on the impulse response

% .. then re-filter each side with this median average impulse response
dinsf = [conv(f1,dins((xx+1):end,1)),conv(f2,dins((xx+1):end,2))];

% .. which we can subtract out
ll = min(length(dmr),length(dinsf));
dvx = dmr(1:ll,:) - dinsf(1:ll,:);
soundsc(dvx(ix,:),sr);

% You can do OK with weiner enhancement even without cancelation
fftlen = 2048;
[mixwf,M] = weinerenhance(dmr, dins, 12.0, 2.0, fftlen);
ff = [0:fftlen/2]*sr/fftlen;
tt = [1:size(M,2)]*fftlen/4/sr;
imagesc(tt,ff,M(:,:,1)); axis xy  % the spectrogram mask
xlabel('Time'); ylabel('Frequency');
axis([20 40 0 4000])
soundsc(mixwf(ix,:), sr);
% but it sounds better based on the enhanced version
[reswf,M] = weinerenhance(resid, targ, 12.0, 2.0);
soundsc(reswf(ix,:), sr);

%% Command line version
%
% remixavier.m wraps these processes into a single function,
% suitable for turning into a compiled Matlab command-line binary:

remixavier -mix ../Data/mc-paul-mix.mp3 -ref ../Data/mc-paul-instr.mp3 -out tmp.wav -dur 60 -weiner_thresh 3.0 -gain 0.9
% use -help to see all the options
remixavier -help


%% Still to do
%
% When the reference signal has very low energy, the coupling
% estimation goes crazy trying to boost it up to get rid of some of
% the energy.  We should put in some kind of
% regularization/threshold to stop this.
%
% We don't expect the coupling filter to vary much along time, so
% we ought to be able to get an improvement by smoothing it along
% time (as in the median filtering on the Imogen Heap example).
% However, if there is any clock drift, we can't assume
% sample-level alignment of the individual impulse response
% estimates.  We could, however, estimate a single timing
% difference between each pair of impulse responses, then average
% them after backing that out.  For instance, we could fit a linear
% phase model to the phase responses of each individual coupling
% IR, then average their zero-phase versions, then reintroduce the
% individual phases (delays) to redistribute over each segment.

%% See also
%
% This project was developed in collaboration with Colin Raffel as part of 
% <http://labrosa.ee.columbia.edu/hamr2013/ HAMR 2013>.
% There is another page describing the project in the 
% <http://labrosa.ee.columbia.edu/hamr2013/proceedings/doku.php/remixavier HAMR Proceedings - Remixavier>
% The code and data are on github:
% <https://github.com/craffel/remixavier>

%% Installation
% 
% This package has been compiled for several targets 
% using the Matlab compiler.  You will also need 
% to download and install the Matlab Compiler Runtime (MCR) Installer. 
% Please see the table below:
%
% <html>
% <table border=1>
% <tr><th>Architecture</th><th>Compiled package</th><th>MCR Installer</th></tr>
% <tr><td>64 bit Linux</td>
% <td><a href="remixavier_GLNXA64.zip">remixavier_GLNXA64.zip</a></td>
% <td><a href="http://www.ee.columbia.edu/~dpwe/tmp/MCRInstaller_glnxa64.bin">Linux 64 bit MCR Installer</a></td></tr>
% <tr><td>64 bit MacOS</td>
% <td><a href="remixavier_MACI64.zip">remixavier_MACI64.zip</a></td>
% <td><a href="http://www.ee.columbia.edu/~dpwe/tmp/MCRInstaller.dmg">MACI64 MCR Installer</a></td></tr>
% </table></html>
% 
% The original Matlab code used to build this compiled target is 
% available at <http://www.ee.columbia.edu/~dpwe/resources/matlab/remixavier>
%
% All sources are in the package <remixavier-v@VER@.zip>.
%
% Feel free to contact me with any problems.

%% Notes
%
% The included function <audioread.m audioread> is able to read a
% wide range of sound file types, but relies on a number of other
% packages and/or support functions being installed.  Most obscure
% of these is  ReadSound, a MEX wrapper I wrote for the dpwelib
% sound file interface.  See the 
% <http://labrosa.ee.columbia.edu/matlab/audioread/ audioread homepage>
% for more details.

%% Changelog

% v0.01  2013-07-01 Initial release

% Last updated: $Date: 2011/12/09 20:30:34 $
% Dan Ellis <dpwe@ee.columbia.edu>
