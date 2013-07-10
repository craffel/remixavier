function [y,a,b] = deskew(dr,dp,sr,deskew_sr,dosquare,PLOT)
% [y,o,k] = deskew(ref,part,sr,deskew_sr,dosquare,PLOT)
%    y is a version of ref that temporally aligns as well as possible
%    with part, which involves delaying the start by o and stretching
%    the timebase by factor k.
%    xcorr performed at sampling rate deskew_sr (default 1000).
%    if dosquare nonzero, square signals before xcorr (default 1).
% 2013-06-29 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3; sr = 44100; end
if nargin < 4; deskew_sr = 0; end
if nargin < 5; dosquare = 1; end
if nargin < 6; PLOT = 0; end

if deskew_sr <= 0;  deskew_sr = 1000;  end


% convert to mono
if size(dr,2) > 1;  drm = sum(dr,2); else   drm = dr; end
if size(dp,2) > 1;  dpm = sum(dp,2); else   dpm = dp; end

% This is a cut down version of skewview
% Set the skewview parameters
xcorrwinsec = 10.0;
xcorrhopsec = 2.0;
xcorrmaxlagsec = 2.0;
backoffsec = 0.5;
xcorrpeakthresh = 0.2;
fitthresh = 2.0;

% Find best global xcorr to nearest 1 ms
%dosquare = 1;
n = find_skew(drm, dpm, [], round(sr/deskew_sr), dosquare);
initialdelay = n/sr;
if n < 0
  msg = ' (part starts before ref)';
elseif n > 0
  msg = ' (ref starts before part)';
else
  msg = '';
end

disp(['Inital estimate of t_ref - t_part = ', ...
      sprintf('%.6f',initialdelay), msg]);

% apply backoff so final time offset is most likely to be a trim,
% not insert
n = n - round(backoffsec*sr);
initialdelay = n/sr;  %%%% reassignment

if n > 0
  % chop off first part of ref
  drm = drm((n+1):end);
else
  % dr starts midway through dp, pre-pad it with zeros to compensate
  drm = [zeros(-n,1);drm];
end

xcorrwin = round(sr * xcorrwinsec);
xcorrmaxlag = round(sr * xcorrmaxlagsec);
xcorrhop = round(sr * xcorrhopsec);

disp('Calculating short-time cross-correlation...');
[Z,E] = stxcorr(drm,dpm,xcorrwin,xcorrhop,xcorrmaxlag);
% normalized xcorr
ZN = Z.*repmat(1./E,size(Z,1),1);

if -min(ZN(:))>max(ZN(:))  % -ve correlation dominates
  zsgn = -1;
else
  zsgn = 1;
end

[zmax,zmaxpos] = max(zsgn*ZN);

% remove points where correlation is much lower than peak
zmaxpos(find(zmax<(xcorrpeakthresh*max(zmax)))) = NaN;

% actual times that corresponds to
zmaxsec = (zmaxpos-xcorrmaxlag-1)/sr;
% hence best linear fit?
tt = [1:length(zmaxpos)]*xcorrhop/sr;
[a,b,lfs,lfp] = linfit(tt, zmaxsec, fitthresh); %,1 for debug in linfit
a = a+1;
disp(sprintf(' Lin fit stats: sd = %.6f prop pts = %.3f',lfs,lfp));
disp(sprintf(' Lin fit: t_ref = %.6f t_part + %.3f', a,b+initialdelay));

if PLOT
  ll = initialdelay + [-xcorrmaxlag:xcorrmaxlag]/sr;

  imagesc(tt,ll,ZN); axis('xy');
  colormap(1-gray);
  colorbar

  xlabel('t_part / sec','interpreter','none')
  ylabel('t_ref - t_part / sec','interpreter','none');
  
  hold on; 
  plot(tt, initialdelay + zmaxsec,'.y'); % yellow dots
  plot(tt, (a-1)*tt + b + initialdelay, '-y');
  hold off;
end

% Apply time offset as pad/trim
n2 = round(b*sr);
if n2 >= 0
  drm = drm((n2+1):end);
else
  drm = [zeros(-n2,1);drm];
end

% total skew
n = n + n2;
% apply to stereo input
if n > 0
  dr = dr(n+1:end,:);
else
  dr = [zeros(-n,size(dr,2));dr];
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
  disp(['Resampling ratio: ',sprintf('%d/%d=%.6f',p(xx),q(xx),p(xx)/q(xx))]);
  for i = 1:size(dr,2)
    y(:,i) = resample(dr(:,i),q(xx),p(xx));
  end
else
  % too close to straight time to resample
  y = dr;
end
