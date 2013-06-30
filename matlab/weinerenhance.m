function y = weinerenhance(target, accomp, thresh, transit, fftlen)
% y = weinerenhance(target, accomp, thresh, transit, fftlen)
%   Enhance a target by keeping only cells where it is large
%   compared to an estimated accompaniment (e.g. the resid and targ
%   from find_in_mix).
%   thresh is the dB threshold (-6.0 by default) and transit is the
%   sigmoid transition width (3.0).  fftlen is the length of the
%   underlying fft (2048).
% 2013-06-30 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3;  thresh = -6.0; end
if nargin < 4;  tranist = 3.0; end
if nargin < 5;  fftlen = 2048; end

DTARG = specgram(target, fftlen);
DCOMP = specgram(accomp, fftlen);
DBRAT = 20*log10(abs(DTARG)) - 20*log10(abs(DCOMP));

% It might be worth temporally smoothing DBRAT here to reduce
% artefacts

% Generate a mask by putting the dB difference through a sigmoid
% with a midpoint at -6 dB, and a transition width ~ 3 dB
y = ispecgram(DTARG .* sigmoid(DBRAT,-6,3));
