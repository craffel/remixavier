function [T,A,H,E,WFs] = decomp_lin_win(Y,X,L,W, doplot)
% [T,A,Hs,E,WFs] = decomp_lin_win(Y,X,L,W)
%     Decompose Y into a component T that is a linear filter
%     applied to X, and a residual A = Y - T.  The estimated linear 
%     filter is L points long (and causal).  Do this repeatedly for 
%     overlapping segments of Y that are W points long, and
%     overlap-add the results.  Hs returns a matrix of filter estimate
% 2010-12-01 Dan Ellis dpwe@ee.columbia.edu for RATS SNR

if nargin < 5; doplot = 0; end

lenY = length(Y);

% Force W to be even
W = 2*round(W/2);

overlapfct = 2;

hop = round(W/overlapfct);

% how many windows?
numwin = 1+ ceil((lenY-W)/hop);

% pad inputs
lastpoint = round(W+(numwin-1)*hop);
X(lastpoint) = 0;
Y(lastpoint) = 0;

% all signals are single columns
T = zeros(lastpoint,1);
A = zeros(lastpoint,1);

hpad = W/2 - hop;
RW = [zeros(hpad,1);hann(2*hop);zeros(hpad,1)];
RWst = RW;
RWen = RW;
RWst(1:(length(RWst)/2)) = 1;
RWen((length(RWst)/2+1):end) = 1;

H = zeros(L,numwin);
E = zeros(1,numwin);

WFs = zeros(2*L+1,numwin);

for i = 1:numwin
  ss = (i-1)*hop+[1:W];
  % record per-window energies
  E(i) = std(X(ss));

  [TT,AA,HH,RXY,WF] = decomp_lin(Y(ss),X(ss),L,doplot);
  if i == 1
    WW = RWst;
  elseif i == numwin
    WW = RWen;
  else
    WW = RW;
  end
  T(ss) = T(ss) + WW.*TT;
  A(ss) = A(ss) + WW.*AA;
  H(:,i) = HH;
  WFs(:,i) = WF';
  
  if doplot
  
    fftlen = 128;
    sr = 8000;
    fmax = 4000;
    subplot(411)
    specgram(Y(ss),fftlen,sr);
    title('mix');
    axis([0 length(ss)/sr 0 fmax]);
    caxis([-40 20]);
    subplot(412)
    specgram(X(ss),fftlen,sr);
    title('clean');
    axis([0 length(ss)/sr 0 fmax]);
    caxis([-40 20]);
    subplot(413)
    specgram(TT,fftlen,sr);
    title('targ');
    axis([0 length(ss)/sr 0 fmax]);
    caxis([-40 20]);
    subplot(414)
    specgram(AA,fftlen,sr);
    title('noise');
    axis([0 length(ss)/sr 0 fmax]);
    caxis([-40 20]);

    pause;
  
  end
  
end

% Trim outputs to be like input
T = T(1:lenY);
A = A(1:lenY);
