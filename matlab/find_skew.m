function [n,xc,ll] = find_skew(orig, part, range)
% [n,xc,ll] = find_skew(orig, part, range)
%    <orig> is an original waveform; <part> is a candidate copy.  
%    <range> (default: length(part)/2) specifies the absolute maximum value
%    to search for <n>, or, if a two-element vector, the min and max values
%    to try for <n>.
%    Return in <n> the optimal alignment within <orig> of the start of
%    <part> based on maximum cross-correlation (can be negative if <part>
%    includes material before beginning of <orig>).
%    <xc> returns the actual cross-correlation function, and <ll>
%    the corresponding lag indices
% 2011-02-11 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3;  range = []; end

if length(range) == 0;
  rangemin = -length(part)+1;
  rangemax = length(orig);
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
%%xc = xcorr(orig, part, maxlag);
%%xcorig = find(abs(xc)==max(abs(xc)))-(maxlag+1);
%xc = xc(rangevals + maxlag + 1);

% cross-correlate by convolution
xcr = fftfilt(flipud(part),[orig;zeros(length(part),1)]);
% actual lags this result corresponds to
%ll = (-length(part))+[1:length(xcr)];
lagmin = -length(part)+1;

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
ll = [rangemin:rangemax]';

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

%disp(['Optimal alignment: part is delayed by ',num2str(xcmax),' samples']);

n = xcmax;
