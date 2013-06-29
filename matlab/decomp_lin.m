function [s_target, e_artif, H, RYX, WF]=decomp_lin(Y,X,L, doplot)
% [s_target, e_artif, H]=decomp_lin(Y,X,L)
%    Y is a corrupted version of X.  Estimate the linear filter H
%    of order L that best predicts X in Y, and return two components
%    s_target = conv(X,H) and e_artif = Y - s_target.
% 2010-12-01 Dan Ellis dpwe@ee.columbia.edu  for RATS SNR analysis

if nargin < 4; doplot = 0; end

% Apply raised cosine ramps at both ends to avoid ringing artefacts
ramptime = 50;
win = 0.5*(1-cos([0:ramptime-1]/ramptime * pi));
wwin = [win';ones(length(Y)-2*ramptime,1);fliplr(win)'];


% Construct whitening filter for clean signal
do_whiten = 1;

  whiteord = 2*L;

  [XW,WF,GF] = whiten(wwin.*X,whiteord);

  % apply gain
  WF = WF*GF;
  
if do_whiten
  % Whiten the output signal
  YW = filter(WF,1,wwin.*Y);
else
  XW = wwin.*X;
  YW = wwin.*Y;
end

% Cross-correlate
RYX = xcorr(YW,XW,L-1);

% Estimate the (causal)linear coupling filter
H = RYX(L-1+[1:L]);

% .. so the linear portion of X in Y would be 
s_target = filter(H,1,X);

% .. and the nonlinear residual is:
e_artif = Y - s_target;

if doplot
  subplot(421); plot(X); title('clean');
  subplot(423); plot(Y); title('noisy');
  subplot(425); plot(e_artif); title('residual');
  subplot(427); plot(-(L-1):(L-1),RYX); title('xcorr'); 
  axis([-200 200 -1 1]); grid

  subplot(422); plot(X); title('clean'); axis([0 500 -.1 .1])
  subplot(424); plot(1:length(Y),Y,1:length(s_target),s_target); ...
      title('noisy + pred'); axis([0 500 -.1 .1])
  subplot(426); plot(e_artif); title('residual'); axis([0 500 -.1 .1])
  subplot(428); plot(WF); title('white filt'); axis([0 500 -5 5])

  pause
end
