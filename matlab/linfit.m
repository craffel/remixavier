function [A,B,S,P] = linfit(X,Y,T,VIZ)
% [a,b,s] = linfit(x,y,t,viz)
%   Try to find a linear regression of y onto x
%   s.t. y[n] = a x[n] + b
%   Do this with rejection of outliers more than t x the 10th
%   percentile (default 2) of distances away
%   from the fit (i.e. where |y - ax + b| > t).  
%   Iterative, heuristic solution.
%   <s> returns the SD of included points around lin fit.
%   <p> returns the proportion of points included in lin fit.
% 2012-09-24 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3;  T = 2.0;  end
if nargin < 4;  VIZ = 0;  end

nx = length(X);

% Initial estimate
nits = 0;
keep = find(~isnan(X) & ~isnan(Y));
nkeep = length(keep)+1;

while nits < 20 && length(keep) < nkeep
  nkeep = length(keep);
  
  % slope as the average first-order difference
%  dXk = diff(X(keep));
%  dYk = diff(Y(keep));
%  kk = find(dXk ~= 0);
%%  A = sum(diff(X(keep)).*(diff(Y(keep))./diff(X(keep))))/sum(diff(X(keep)));
%%  A = sum(dXk(kk).*(dYk(kk)./dXk(kk)))/sum(dXk(kk));
%  % mode?
%%  [NN,XX] = hist(dYk(kk)./dXk(kk),length(kk)/3);
%%  [vv,xx] = max(NN);
%%  A = XX(xx);
%  % median
%  A = percentile(dYk(kk)./dXk(kk),0.5);
%  % offset as the residual
%%  B = mean(Y(keep) - A.*X(keep));
%  B = percentile(Y(keep) - A * X(keep),0.5);

  % Least squares
  C = [X(keep)',ones(length(keep),1)];
  AB = inv(C'*C)*C'*(Y(keep)');
  A = AB(1);
  B = AB(2);


  % Refine
  D = Y - (A * X + B);

  Dth = max(T * percentile(abs(D),0.2), 0.001); % don't be tighter than 1ms

  oldkeep = keep;
  keep = find(abs(D) < Dth);

  % Also remove points whose slopes are extreme
  slopes = diff(X(keep))./diff(Y(keep));
  [vv,xx] = sort(slopes);
  % remove first and last deciles
  decile = round(0.1*length(slopes));
  keep = keep(xx((decile+1):(end-decile)));
  % Will remove the point immediately preceding an extreme slope,
  % but outliers will have large slopes both before and after, so
  % will tend to be taken.  
  
  % Proportion of points retained
  P = length(keep)/length(X);
  % SD around line of retained points
  S = sqrt(mean( (Y(keep)-(A*X(keep)+B)).^2));
  
  if VIZ
    disp(['Keep from ',num2str(length(oldkeep)),' to ', ...
          num2str(length(keep)), ' points']);
    plot(X,Y,'.c');
    hold on
    plot(X(oldkeep),Y(oldkeep),'.b');
    plot(X,A*X+B,'-r');
    plot(X(keep),Y(keep),'.r');
    hold off
    pause
  end
  
  nits = nits + 1;
end
