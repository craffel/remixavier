function [wavout,sr] = remixavier(varargin)
% [wavout,sr] = remixavier(varargin)
%   Utility to separate and enhance certain sources when provided
%   with alternate mixes that do and do not include them.
%
% Usage;
%   remixavier [options]
%
%   See http://labrosa.ee.columbia.edu/projects/remixavier
%   for full documentation
%
% 2013-07-01 Dan Ellis dpwe@ee.columbia.edu
% $Header: $

VERSION = 0.03;
DATE = 20130709;
progname = ['*** remixavier v', num2str(VERSION),' of ', num2str(DATE)];

[P,X] = get_params(varargin, ...
               {'mix', '', 'audio with extra source(s)', ...
                'part', '', 'audio without extra source(s)', ...
                'out', '', 'write audio output to this file', ...
                'alignout', '', 'write aligned part to this file', ...
                'gain', 1, 'scale output by this to avoid clipping', ...
                'mono', 0, 'force files to mono?', ...
                'flip_stereo', 0, 'flip L and R of part', ...
                'samplerate', 0, 'resample inputs to this rate', ...
                'mix_start', 0, 'start reading mix file from this point', ...
                'part_start', 0, 'start reading part file from this point', ...
                'dur', 0, 'limit processing to this duration', ...
                'hpf_frq', 0, 'high-pass frequency break frequency', ...
                'ir_dur', 0.015, 'total time extent of coupling filter', ...
                'ir_pre', 0.005, 'pre-echo time in coupling filter', ...
                't_win', 1, 'duration of filter estimation window in sec', ...
                't_hop', 0.5, 'hop between successive estimation wins in sec', ...
                'deskew_its', 1, 'how many times to pass through deskew', ...
                'deskew_sr', 1000, 'sampling rate for initial deskew xcorr', ...
                'square', 1, 'square waveforms before initial xcorr', ...
                'do_plot', 0, 'plot the results of deskewing', ...
                'wiener_win', 0.050, 'STFT duration for Wiener filter', ...
                'wiener_thresh', -Inf, 'local SNR threshold (in dB) for Wiener enhancement', ...
                'wiener_width', 3.0, 'transition width (in dB) for Wiener'}, ...
               progname);

if isempty(P)
  % error/help from get_params
  if nargout > 0
    wavout = []; sr = []; 
  end
  return;  % exit the program
end


if length(P.mix) == 0 || length(P.part) == 0
  error('remixavier: you must specify both -mix and -part');
end

% Read in both files
[dmix, sr] = audioread(P.mix, P.samplerate, P.mono, P.mix_start, P.dur);
[dpart, sr] = audioread(P.part, sr, P.mono, P.part_start, P.dur);

dmix = P.gain * dmix;

if P.hpf_frq > 0
  % design filter
  order = 8;
  ripple = 2.0; % dB
  Wp = P.hpf_frq / (sr/2);
  [b,a] = cheby1(order,ripple,Wp,'high');
  dmix = filter(b,a,dmix);
  dpart = filter(b,a,dpart);
end

if (P.flip_stereo == 1) && (size(dpart, 2) == 2)
  disp('stereo flip')
  dpart = dpart(:,[2 1]);
end

% Apply deskewing (time offset and clock drift correction)
% (May be run multiple times for better results)
do_plot = P.do_plot;
for i = 1:P.deskew_its
  dpart = deskew(dpart, dmix, sr, P.deskew_sr, P.square, do_plot);
  do_plot = 0;  % only plot on first iteration
end

if length(P.alignout) > 0
  % Maybe write out aligned version
  audiowrite(dpart, sr, P.alignout);
  disp(['Aligned part written to ',P.alignout]);
end

% Estimate coupling filters & cancel for each channel individually
for i = 1:size(dmix,2)
  partc = 1+mod(i-1, size(dpart,2));  % cycle part chans if too few
%  [resid(:,i), targ(:,i), filt, SNR, del, filts] = ...
  [res, tar, filt, SNR, del, filts] = ...
      find_in_mix(dmix(:,i), dpart(:,partc), sr, ...
                  P.ir_dur, P.ir_pre, ...
                  P.t_win, P.t_hop);
  if i == 1
    resid = res; targ = tar;
    lr = length(resid);
  else
    res(lr+1) = 0;  % pad with zeros if too short
    tar(lr+1) = 0;
    resid(:,i) = res(1:lr); targ(:,i) = tar(1:lr);  % take first
                                                    % bit if too long
  end
end

if P.do_plot > 1
  % Plot the time-local coupling filters (last channel only)
  % filter IR time base
  tt = [1:size(filts,1)]/sr;
  % times of individual short-time window (every 4 sec)
  tw = P.t_hop*[1:size(filts,2)];
  % plot
  imagesc(tt,tw,filts'); axis('xy');
  xlabel('time / sec')
  ylabel('window time / sec')
  title('local coupling filter impulse responses (part -> mix)')
  % scale down impulse response extremes
  caxis([-1 1]*min(2, max(abs(caxis()))));
  colorbar;
end

% Apply Wiener enhancement
if P.wiener_thresh > -Inf
  fftlen = 2^round(log(sr*P.wiener_win)/log(2));  % nearestpow2
  resid = wienerenhance(resid, targ, P.wiener_thresh, P.wiener_width, ...
                        fftlen);
end

% Save result
if length(P.out) > 0
  audiowrite(resid, sr, P.out);
  disp(['Canceled audio written to ',P.out]);
end
% Maybe return result
if nargout > 0
  wavout = resid;
end

if P.do_plot == 3
  tfft = 0.032;
  fftlen = 2^round(log(sr*tfft)/log(2));
  p1 = subplot(311);
  specgram(sum(dmix,2), fftlen, sr);
  cax = max(caxis());
  title(['mix - ',P.mix], 'interpreter', 'none');
  p2 = subplot(312);
  specgram(sum(targ,2), fftlen, sr);
  cax = max([cax, caxis()]);
  title(['part - ',P.part], 'interpreter', 'none');
  p3 = subplot(313);
  specgram(sum(resid,2), fftlen, sr);
  cax = max([cax, caxis()]);
  title(['residual'], 'interpreter', 'none');
  subplot(311)
  caxis(cax+[-80 0]);
  subplot(312)
  caxis(cax+[-80 0]);
  subplot(313)
  caxis(cax+[-80 0]);
  linkaxes
  % figure gated SNR
  Emix = sum(abs(specgram(sum(dmix,2), fftlen, sr)).^2);
  Etrg = sum(abs(specgram(sum(targ,2), fftlen, sr)).^2);
  Ersd = sum(abs(specgram(sum(resid,2), fftlen, sr)).^2);
  
  lowthDB = -30;
  lowth = 10^(lowthDB/20);
  goodcols = find(Etrg > lowth*max(Etrg)); 
  segSNR = 10*log10(mean(Etrg(goodcols)./Ersd(goodcols)));
  disp(['segSNR = ', num2str(segSNR),' (over ', ...
        num2str(length(goodcols)/length(Etrg)*100),'% of frames)']);  

  subplot(312); hold on; plot(goodcols*(fftlen/2)/sr, 200, '.m'); hold off
end
