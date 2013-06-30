% demo of matlab version of remixavier
% 2013-06-28 Dan Ellis dpwe@ee.columbia,edu

% Load in mix and acapella as mono files
% These tracks diverge at the end, so just work on the first minute
sr = 44100;
[dmix,sr] = mp3read('../Data/mc-paul-mix.mp3',[0 60*sr],1);
[dins,sr] = mp3read('../Data/mc-paul-instr.mp3',[0 60*sr],1);

% Attempt to trim and resample the full version to line up as well
% as possible with the acapella
dmr = deskew(dmix, dins);
% It gets better when you repeat it
dmr = deskew(dmr, dins);
% resampling can't handle ratios below 10 ppm, will just skip
% beyond that.

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
% times of individual short-time window
tw = 0.25*[1:size(filts,2)];
% plot
imagesc(tt,tw,filts'); axis('xy');
xlabel('time / sec')
ylabel('window time / sec')
title('local coupling filter impulse responses (cap -> mix)')
% scale down impulse response extremes
caxis([-2 2])


% The duffy track has the vocals in stereo, we can cancel left and
% right separately to good effect
[dmix,sr] = mp3read('../Data/Duffy.WarwickAvenue.mp3');
[dcap,sr] = mp3read('../Data/duffy_-_warwick_avenue_acapella.mp3');

% Deskew will process stereo files.  Skew is estimated from an
% internally-generated mono mix
dmr = deskew(dmix, dcap);
dmr = deskew(dmr, dcap);

clear resid

for i = 1:size(dmr,2)
  tic; [resid(:,i), targ, filt, SNR, del, filts] = ...
                find_in_mix(dmr(:,i),dcap(:,i),sr,0.006,0.003); toc
end
  
soundsc(resid(ix,:), sr);


% The Imogen Heap track has no time skew, so we can average
% together all the local filter responses - or median is better
[dmix,sr] = audioread('../Data/10-Aha_.m4a');
[dins,sr] = audioread('../Data/23-Aha_Instrumental_Version_.m4a');

% Deskew once just to line them up in time
dmr = deskew(dmix, dins);

clear resid filts

% Align each channel, and store all the individual filters
for i = 1:size(dmr,2)
  tic; [resid(:,i), targ, filt, SNR, del, filts{i}] = ...
                find_in_mix(dmr(:,i),dins(:,i),sr,0.006,0.003); toc
end

soundsc(resid(ix,:), sr);

% but form a grand average filter for each side
f1 = median(filts{1}');
f2 = median(filts{2}');

% It has a pre-echo, so trim that from the convolution
[vv,xx] = max(abs(f1));
% xx is the index of the peak on the impulse response

% .. then re-filter each side with this median average impulse response
dinsf = [conv(f1,dins((xx+1):end,1)),conv(f2,dins((xx+1):end,2))];

% .. which we can subtract out
ll = min(length(dmr),length(dinsf));
dvx = dmr(1:ll,:) - dinsf(1:ll,:);
soundsc(dvx(ix,:),sr);


% Message In A Bottle is an ideal case - plain subtraction of mix
% and instrumental yeilds clean vocals.  But how do we do?

sr = 44100;
dmix = mean(wavread('../Data/message-in-a-bottle-mix.wav'),2);
dins = mean(wavread('../Data/message-in-a-bottle-ins.wav'),2);

dmr = deskew(dmix, dins);

tic; [resid, targ, filt, SNR, del, filts] = ...
      find_in_mix(dmr,dins,sr,0.013,0.003); toc

soundsc(resid(ix,:), sr);

% We can apply a "wiener filter" (scaling of spectrogram magnitude
% cells) to further reduce residual artifacts.  In particular, we
% can suppress cells where the energy in the estimated vocals is
% significantly lower than the energy in the instrumental line
% projected into the mix.  weinerenhace with a threshold so that
% energy in the residual that is below -6 dB of the accompaniment
% is suppressed

reswf = weinerenhance(resid, targ, -6.0);
soundsc(reswf(ix,:), sr);

% We can measure SNR by canceling against the true vocals, which
% are simply the difference of dmix and dins (for this perfect example)
dvox = dmix - dins;
soundsc(dvox(ix,:), sr);  % Yes, sounds clean
[r2, t2, f2, S2, d2, fs2] = find_in_mix(resid,dvox,sr,0.010,0.003); 
%Delay = 0.000000 s
%SNR = 19.9197 dB
[r2, t2, f2, S2, d2, fs2] = find_in_mix(reswf,dvox,sr,0.010,0.003);
%Delay = 0.000000 s
%SNR = 16.1096 dB
% Weiner filtering introduces more artefact energy than it removes
% interference.
