function [y,a,b] = deskew(dm,dr,sr)
% [y,o,k] = deskew(x,r,sr)
%    y is a version of x that temporally aligns as well as possible
%    with r, which involves delaying the start by o and stretching
%    the timebase by factor k.
% 2013-06-29 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3; sr = 44100; end

PLOT = 0;

% Cut down version of skewview
xcorrwinsec = 4.0;
xcorrhopsec = 1.0;
xcorrmaxlagsec = 2.0;
backoffsec = 0.5;
xcorrpeakthresh = 0.1;
fitthresh = 2.0;

% Find best global xcorr
[n,xc] = find_skew(dr, dm);
initialdelay = n/sr;
disp(['Initial delay = ',num2str(initialdelay)]);

% apply backoff so final time offset is most likely to be a trim,
% not insert
n = n - round(backoffsec*sr);
initialdelay = n/sr;  %%%% reassignment

if n > 0
  % dm starts midway through dr, chop off first part of dr
  %dr = dr((n+1):end);
  dm = [zeros(n,1);dm];
else
  % chop off first part of mix
  dm = dm(((-n)+1):end);
end

xcorrwin = round(sr * xcorrwinsec);
xcorrmaxlag = round(sr * xcorrmaxlagsec);
xcorrhop = round(sr * xcorrhopsec);

disp('Calculating short-time cross-correlation...');
[Z,E] = stxcorr(dr,dm,xcorrwin,xcorrhop,xcorrmaxlag);
% normalized xcorr
ZN = Z.*repmat(1./E,size(Z,1),1);

[zmax,zmaxpos] = max(ZN);

% remove points where correlation is much lower than peak
zmaxpos(find(zmax<(xcorrpeakthresh*max(zmax)))) = NaN;

% actual times that corresponds to
zmaxsec = (zmaxpos-xcorrmaxlag-1)/sr;
% hence best linear fit?
tt = [1:length(zmaxpos)]*xcorrhop/sr;
[a,b] = linfit(tt, zmaxsec, fitthresh); %,1 for debug in linfit
a = a+1;
disp(sprintf('Lin fit: y = %.6f x + %.3f', a,b+initialdelay));

if PLOT
  ll = initialdelay + [-xcorrmaxlag:xcorrmaxlag]/sr;

  imagesc(tt,ll,ZN); axis('xy');
  colormap(1-gray);
  colorbar

  hold on; 
  plot(tt, initialdelay + zmaxsec,'.r');
  hold off;
end

% Apply time offset as pad/trim
if b < 0
  dm = dm((round(b*sr)+1):end);
else
  dm = [zeros(round(b*sr),1);dm];
end

% Apply time scaling via resampling
% Find p/q s.t. p/q ~= a
% approx 1 part in 1/(a-1)
p0 = floor(1/abs(a-1));
if p0 < 2^15
  % exhaustive search for pair of integers closest to desired rate
  p = p0:(2^15);
  q = round(p./a); 
  er = (a - p./q); 
  [ee,xx] = min(abs(er));
  disp(sprintf('%d/%d=%.6f',p(xx),q(xx),p(xx)/q(xx)));
  y = resample(dm,p(xx),q(xx));
else
  % too close to straight time to resample
  y = dm;
end
