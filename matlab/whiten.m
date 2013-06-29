function [Y,A,G] = whiten(X,P)
% y = whiten(x,p)
%  Fit a p'th order LPC to the whole of X, and inverse-filter by
%  it.
% 2010-11-27 Dan Ellis dpwe@ee.columbia.edu

A = lpc(X,P);
% handle degenerate case if X is all zeros
A(isnan(A)) = 0;
Y = filter(A,1,X);
if norm(Y) > 0
  G = 1/norm(Y);
else
  G = 1;
end

Y = G*Y;

