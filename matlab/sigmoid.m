function y = sigmoid(x,t,s)
% y = sigmoid(x,t,s)
%   Simple sigmoid function goes from 0 to 1 as x goes from t-s to
%   t+s.
% 2013-06-29 Dan Ellis dpwe@ee.columbia.edu

xd = (x-t)/s;

y = 0.5+0.5*(xd ./ sqrt(1 + xd.^2));