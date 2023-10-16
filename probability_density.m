function [pd,F,n,Fdot_mean] = probability_density(tbl,dF,plotting)
% Histogram of probability densities from experimennts
% Input:
%  tbl:  Tu or Tr from unfolding analysis
%  dF:   Force value spacing
%  plotting: 0 ( default) NO plotting
%            1 to display the histogram in new plot
% Output: 
%  pd: Column vector of probality desities.  The probability of unfolding 
%      (or refolding) in the interval <F(i)-dF,F(i+dF)< is pd(i)*dF
%  n: nUmber of rows in tbl
%  F: Force vector (useful if F is not specified in the call)
%  Fdot_mean: mean value of Fdot (= df/dt, where f(t) is the pulling force)

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
  Fdot_mean = mean(tbl.Fdot);  
  n = size(tbl,1);
  edges = (F(1)-0.5*dF):dF:(F(end)+0.5*dF);
  Values = histcounts(tbl.Force,edges);
  pd = Values'/sum(Values)/dF;  % Probability density
  if plotting
    figure;
    histogram(tbl.Force,edges);
    cla;
    b = bar(F,pd,1); 
    set(b,'facecolor',[0.2,.7,.9])
    title(sprintf('%d events, mean temperaure = %5.1f°C',n,mean(tbl.Temperature)))
    ylabel('Probability density')
    xlabel('Force (pN)')
  end
end

