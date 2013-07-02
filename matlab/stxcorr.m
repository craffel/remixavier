function [Z,E] = stxcorr(X,Y,W,H,N,OLDWAY)
% [Z,E] = stxcorr(X,Y,W,H,N)  short-time cross-correlation
%    X and Y are segmented into W-point stretches, hopped by H samples. 
%    The cross-correlation of these, out to +/- N points, are returned 
%    as columns of Z.
%    Result is timing of X relative to Y, so that if X is equal to
%    some padding samples followed by Y (i.e., Y delayed) then the
%    cross correlation has peak values beyond its midpoint (in each
%    column of Z).
%    E returns theoretical maximum (harmonic mean energy) for each
%    col of Z.
% 2004-04-28 dpwe@ee.columbia.edu

if nargin < 2;   Y = []; end
if nargin < 3;   W = 2000;  end
if nargin < 4;   H = round(W/2); end
if nargin < 5;   N = W; end
if nargin < 6;   OLDWAY = 0; end  % was ismac()

if N > W; N = W; end

npts = 2*N+1;
LX = length(X);
if ~isempty(Y)
  LX = min(LX,length(Y));
end
nfr = 1+floor((LX-W)/H);

Z = zeros(npts, nfr);
E = zeros(1, nfr);

% Hann window
%ww = hann(W); % is not used here!

% time domain version - for a 240k pt signal with W=160 H=N=80,
% took 1.14 sec
if OLDWAY
  for i = 1:nfr
    pts = (i-1)*H + [1:W];
%    WX = W.*X(pts);
%    WY = W.*Y(pts);
%   clearly I was thinking that W was the window, but it isn't
    WX = X(pts);
    WY = Y(pts);
    Z(:,i) = xcorr(WX, WY, N);
    E(i) = sqrt(sum(WX.^2)*sum(WY.^2));
  end
else

  % faster (ahem) version, identical result
  % took 0.16 sec
  % 2012-07-10 Dan Ellis dpwe@ee.columbia.edu
  XX = fft(frame(X(1:LX),W,H),2*W);
  if isempty(Y)
    YY = XX;
  else
    YY = fft(frame(Y(1:LX),W,H),2*W);
  end
  xx = ifft(XX.*conj(YY));
  % This assumes autocorrelation
  %Z = [xx([N+1:-1:2],:);xx([1:N+1],:)];
  %E = Z(N+1,:);

  Z = [xx([end-(N-1):end],:);xx([1:N+1],:)];
  E = sqrt(sum(abs(XX).^2).*sum(abs(YY).^2))/(2*W);
end
