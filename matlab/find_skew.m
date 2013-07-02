function [n,xc,ll] = find_skew(test, ref, range, resolution, dosquare)
% [n,xc,ll] = find_skew(test, ref, range, resolution, dosquare)
%    <test> is a test waveform; <ref> is a reference we are
%    comparing it to.  <n> returns the index of <test> within
%    <ref>; if <test> is ref but with some samples prepended (i.e.,
%    <test> is a delayed version of <ref>), <n> returns as a
%    positive value, the location of the maximum of the
%    cross-correlation function, or equivalently the sample index
%    within <ref> that best lines up with the first sample of <test>.
%    <range> (default: length(ref)/2) specifies the absolute maximum value
%    to search for <n>, or, if a two-element vector, the min and max values
%    to try for <n>.
%    <resolution> is the required accuracy in samples.  Waveforms
%    will be downsampled to do no better than this (0 = no
%    downsampling, default = 16).
%    <dosquare> set to 1 causes cross-correlation to be performed
%    on squared signals (gives better results for signals with
%    small amounts of clock drift).
%    <xc> returns the actual cross-correlation function, and <ll>
%    the corresponding lag indices
% 2011-02-11 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3;  range = []; end
if nargin < 4;  resolution = 16; end
if nargin < 5;  dosquare = 0; end

if resolution > 1
  test = resample(test, 1, resolution);
  ref = resample(ref, 1, resolution);
  range = round(range / resolution);
else
  resolution = 1;
end

% I find I can handle clock drift better if I square the signals
if dosquare
  test = test.^2;
  ref = ref.^2;
end

if length(range) == 0;
  rangemin = -length(ref)+1;
  rangemax = length(test);
elseif length(range) == 1
  rangemin = -abs(range);
  rangemax = abs(range);
else
  rangemin = range(1);
  rangemax = range(2);
end

xc = zeros(rangemax-rangemin+1,1);

%% find offset between clean and noisy
%maxlag = max(abs(range));
%
%%xc = xcorr(test, ref, maxlag);
%%xctest = find(abs(xc)==max(abs(xc)))-(maxlag+1);
%xc = xc(rangevals + maxlag + 1);

% cross-correlate by convolution
xcr = fftfilt(flipud(ref),[test;zeros(length(ref),1)]);
% actual lags this result corresponds to
%ll = (-length(ref))+[1:length(xcr)];
lagmin = -length(ref)+1;

xcmax = min(lagmin + (find(abs(xcr)==max(abs(xcr))))) - 1;

% where in ll does the first value to be returned in xcvals occur?
offset = rangemin - lagmin;   % ll(offset+1) = rangevals(1) or
                              % rangevals(-offset+1) = ll(1)
% which points in ll will we copy?
touse = (max(0,-offset)+1):min(length(xc),length(xcr)-offset);

%xc = xcr;
%ll = range;

xc(touse) = xcr(offset+touse);
%ll = ll(offset+touse);
ll = resolution * [rangemin:rangemax]';

% % trim files to be aligned
% if xcmax > 0
%   dn = dn(xcmax:end);
%   TS = TS + xcmax/sr;
% elseif xcmax < 0
%   dc = dc((-xcmax):end);
% end
% % make files same length
% dlen = min(length(dn),length(dc));
% dn = dn(1:dlen);
% dc = dc(1:dlen);
% TE = TS + dlen/sr;

%disp(['Optimal alignment: ref is delayed by ',num2str(xcmax),' samples']);

n = resolution * xcmax;
