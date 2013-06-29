function [n,xc] = find_skew(orig, part, range)
% [n,xc] = find_skew(orig, part, range)
%    <orig> is an original waveform; <part> is a candidate copy.  
%    Return in <n> the optimal alignment within <orig> of the start of
%    <part> based on maximum cross-correlation (can be negative if <part>
%    includes material before beginning of <orig>).
%    <range> (default: length(part)/2) specifies the absolute maximum value
%    to search for <n>, or, if a two-element vector, the min and max values
%    to try for <n>.  <xc> returns the actual cross-correlation function.
% 2011-02-11 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3;  range = length(part)/2; end

if length(range) == 1
  range = [-range range];
end
rangevals = range(1):range(2);

% find offset between clean and noisy
maxlag = max(abs(range));
xc = xcorr(part, orig, maxlag);
%xcorig = find(abs(xc)==max(abs(xc)))-(maxlag+1);
xc = xc(rangevals + maxlag + 1);
xcmax = rangevals(find(abs(xc)==max(abs(xc))));

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

% disp(['Optimal alignment: part is delayed by ',num2str(xcmax),' samples']);

n = xcmax;