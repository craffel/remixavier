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

VERSION = 0.01;
DATE = 20130701;
progname = ['*** remixavier v', num2str(VERSION),' of ', num2str(DATE)];

[P,X] = get_params(varargin, ...
               {'mix', '', 'audio with extra source(s)', ...
                'ref', '', 'audio without extra source(s)', ...
                'out', '', 'write audio output to this file', ...
                'gain', 1, 'scale output by this to avoid clipping', ...
                'mono', 0, 'force files to mono?', ...
                'samplerate', 0, 'resample inputs to this rate', ...
                'mix_start', 0, 'start reading mix file from this point', ...
                'ref_start', 0, 'start reading ref file from this point', ...
                'dur', 0, 'limit processing to this duration', ...
                'ir_dur', 0.015, 'total time extent of coupling filter', ...
                'ir_pre', 0.005, 'pre-echo time in coupling filter', ...
                't_win', 8, 'duration of filter estimation window in sec', ...
                't_hop', 4, 'hop between successive estimation wins in sec', ...
                'deskew_its', 1, 'how many times to pass through deskew', ...
                'do_plot', 0, 'plot the results of deskewing', ...
                'weiner_win', 0.050, 'STFT duration for Weiner filter', ...
                'weiner_thresh', -Inf, 'local SNR threshold (in dB) for Weiner enhancement', ...
                'weiner_width', 3.0, 'transition width (in dB) for Weiner'}, ...
               progname);

if length(P.mix) == 0 || length(P.ref) == 0
  error('remixavier: you must specify both -mix and -ref');
end

% Read in both files
[dmix, sr] = audioread(P.mix, P.samplerate, P.mono, P.mix_start, P.dur);
[dref, sr] = audioread(P.ref, P.samplerate, P.mono, P.ref_start, P.dur);

% Apply deskewing (time offset and clock drift correction)
% (May be run multiple times for better results)
do_plot = P.do_plot;
for i = 1:P.deskew_its
  dmix = deskew(dmix, dref, sr, do_plot);
  do_plot = 0;  % only plot on first iteration
end

% Estimate coupling filters & cancel for each channel individually
for i = 1:size(dmix,2)
  refc = 1+mod(i-1, size(dref,2));  % cycle ref chans if too few
  [resid(:,i), targ(:,i)] = find_in_mix(dmix(:,1), dref(:,refc), sr, ...
                                        P.ir_dur, P.ir_pre, ...
                                        P.t_win, P.t_hop);
end

% Apply Weiner enhancement
if P.weiner_thresh > -Inf
  fftlen = 2^round(log(sr*P.weiner_win)/log(2));  % nearestpow2
  resid = weinerenhance(resid, targ, P.weiner_thresh, P.weiner_width, ...
                        fftlen);
end

% Save result
if length(P.out) > 0
  audiowrite(P.gain * resid, sr, P.out);
end
% Maybe return result
if nargout > 0
  wavout = resid;
end
