function [pd,F,n,Tmean,Fdot] = probability_density(tbl,dF)
% Histogram of probability densities from experimennts
% Input:
%  tbl:  Tu or Tr from unfolding analysis
%  dF:   Force value spacing
% Output: 
%  pd: Column vector of probality desities.  The probability of unfolding 
%      (or refolding) in the interval <F(i)-dF/2,F(i+dF/2)< is pd(i)*dF
%  F: Force vector 
%  n: nUmber of rows in tbl
%  Tmean: Mean of tbl.Temperature
%  Fdot: mean value of df/dt before unfolding or refolding, 
%     where f(t) is the pulling force)

% Version 1.1: Fdot is now negative for refolding traces

  if nargin < 3
    plotting = 0;
  end
  if nargin < 2
    dF = 1;
  end
  if isempty(dF)
      dF = 1;
  end
  F = (floor(min(tbl.Force)):dF:ceil(max(tbl.Force)))';
  if max(tbl.Pullingspeed)/min(tbl.Pullingspeed) > 4
    warning('Events with a very wide range of pullingspeeds should be analysed separately');
  end
  Fdot = mean(tbl.Fdot);  
  Tmean = mean(tbl.Temperature);
  n = size(tbl,1);
  edges = (F(1)-0.5*dF):dF:(F(end)+0.5*dF);
  Values = histcounts(tbl.Force,edges);
  pd = Values'/sum(Values)/dF;  % Probability density
end

