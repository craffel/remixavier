function [y,a,b] = deskew(dm,dr,sr)
% [y,o,k] = deskew(x,r,sr)
%    y is a version of x that temporally aligns as well as possible
%    with r, which involves delaying the start by o and stretching
%    the timebase by factor k.
% 2013-06-29 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3; sr = 44100; end

% convert to mono
if size(dm,2) > 1;  dmm = sum(dm,2); else   dmm = dm; end
if size(dr,2) > 1;  drm = sum(dr,2); else   drm = dr; end

PLOT = 0;

% Cut down version of skewview
xcorrwinsec = 4.0;
xcorrhopsec = 1.0;
xcorrmaxlagsec = 2.0;
backoffsec = 0.5;
xcorrpeakthresh = 0.1;
fitthresh = 2.0;

% Find best global xcorr
[n,xc] = find_skew(dmm, drm);
initialdelay = n/sr;
disp(['Initial delay = ',sprintf('%.6f',initialdelay)]);

% apply backoff so final time offset is most likely to be a trim,
% not insert
n = n + round(backoffsec*sr);
initialdelay = n/sr;  %%%% reassignment

if n > 0
  % chop off first part of mix
  dmm = dmm((n+1):end);
else
  % dm starts midway through dr, chop off first part of dr
  %dr = dr((-n+1):end);
  dmm = [zeros(-n,1);dmm];
end

xcorrwin = round(sr * xcorrwinsec);
xcorrmaxlag = round(sr * xcorrmaxlagsec);
xcorrhop = round(sr * xcorrhopsec);

disp('Calculating short-time cross-correlation...');
[Z,E] = stxcorr(drm,dmm,xcorrwin,xcorrhop,xcorrmaxlag);
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
disp(sprintf('Lin fit: y = %.6f x + %.3f', a,b-initialdelay));

if PLOT
  ll = -initialdelay + [-xcorrmaxlag:xcorrmaxlag]/sr;

  imagesc(tt,ll,ZN); axis('xy');
  colormap(1-gray);
  colorbar

  hold on; 
  plot(tt, -initialdelay + zmaxsec,'.r');
  hold off;
end

% Apply time offset as pad/trim
n2 = round(-b*sr);
if n2 >= 0
  dmm = dmm((n2+1):end);
else
  dmm = [zeros(-n2,1);dmm];
end

% total skew
n = n + n2;
% apply to stereo input
if n > 0
  dm = dm(n+1:end,:);
else
  dm = [zeros(-n,size(dm,2));dm];
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
  for i = 1:size(dm,2)
    y(:,i) = resample(dm(:,i) ,p(xx), q(xx));
  end
else
  % too close to straight time to resample
  y = dm;
end
