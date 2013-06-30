function [Y,winix] = frame(X,W,H)
% Y = frame(X,W,H)
%   Return Y as a set of columns as W-point segments of X stepped
%   by H. 
%   There's no windowing (i.e. scaling by a tapered window) in here. 
% 2010-11-14 Dan Ellis dpwe@ee.columbia.edu

%lx = length(X);
%nh = 1+ceil((length(X)-W)/H);
%% Pad X to an integral number of windows
%Xp = [X(:)',zeros(1, (W+(nh-1)*H)-lx)];
% No, truncate to have only whole windows
nh = 1+floor((length(X)-W)/H);
Xp = [X(1:W+(nh-1)*H)'];

% Index-fu:
% We build a matrix of indices that pull out the values we want...
% (columns of 1:W, with successive multiples of H added on)
winix = repmat(H*[0:(nh-1)],W,1)+repmat([1:W]',1,nh);

% .. then the output is just the input matrix indexed by these indices
Y = Xp(winix);
