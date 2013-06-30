function [y,M,DTARG] = weinerenhance(target, accomp, thresh, transit, fftlen)
% [y,M,D] = weinerenhance(target, accomp, thresh, transit, fftlen)
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
    [y(:,i),M(:,:,i),DTARG(:,:,i)] = weinerenhance(target(:,i),accomp(:,i), ...
                                                   thresh, transit, ...
                                                   fftlen);
  end
else

  % default, mono case
  DTARG = specgram(target, fftlen);
  DCOMP = specgram(accomp, fftlen);
  frames = min(size(DTARG,2), size(DCOMP,2));
  DBRAT = 20*log10(abs(DTARG(:,1:frames))) - 20*log10(abs(DCOMP(:,1:frames)));

  % It might be worth temporally smoothing DBRAT here to reduce
  % artefacts

  % Generate a mask by putting the dB difference through a sigmoid
  % with a midpoint at -6 dB, and a transition width ~ 3 dB
  M = sigmoid(DBRAT,thresh,transit);
  y = ispecgram(DTARG .* M);

end
