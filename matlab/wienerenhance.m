function [y,M,DTARG] = wienerenhance(target, accomp, thresh, transit, fftlen)
% [y,M,D] = wienerenhance(target, accomp, thresh, transit, fftlen)
%   Enhance a target by keeping only cells where it is large
%   compared to an estimated accompaniment (e.g. the resid and targ
%   from find_in_mix).
%   thresh is the dB threshold (-6.0 by default) and transit is the
%   sigmoid transition width (3.0).  fftlen is the length of the
%   underlying fft (2048).
%   M is the resulting mask, D is the unmodified corresponding STFT
%   of target.
% 2013-06-30 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3;  thresh = -6.0; end
if nargin < 4;  transit = 3.0; end
if nargin < 5;  fftlen = 2048; end

if size(target,2) > 1
  % stereo signal, factor into channels
  for i = 1:size(target,2)
    [y(:,i),M(:,:,i),DTARG(:,:,i)] = wienerenhance(target(:,i),accomp(:,i), ...
                                                   thresh, transit, ...
                                                   fftlen);
  end
else

  % default, mono case
  DTARG = stft(target, fftlen, fftlen, fftlen/4);
  DCOMP = stft(accomp, fftlen, fftlen, fftlen/4);
  frames = min(size(DTARG,2), size(DCOMP,2));
  % remove zero frames (which mess up the log)
  DTARG(find(DTARG==0)) = 1e-5;
  DCOMP(find(DCOMP==0)) = 1e-5;
  % Figure the DB ratio
  DBRAT = 20*log10(abs(DTARG(:,1:frames))) - 20*log10(abs(DCOMP(:,1:frames)));

  % It might be worth temporally smoothing DBRAT here to reduce
  % artefacts

  % Generate a mask by putting the dB difference through a sigmoid
  % with a midpoint at -6 dB, and a transition width ~ 3 dB
  M = sigmoid(DBRAT,thresh,transit);
  y = istft(DTARG(:,1:frames) .* M, fftlen, fftlen, fftlen/4);
  % istft returns as a row
  y = y';
end
