function [pd,F,n] = probdens(force,dF,plotting)
% Histogram of probability densities from experimennts
% About the same as probability_density, but uses force vector as input
%   rather that Matlab Table Tun or Tre
% Input:
%  force:  Array of recorded transition forces
%  dF:   Force value spacing
%  plotting: 0 ( default) NO plotting
%            1 to display the histogram in new plot
% Output: 
%  pd: Column vector of probality desities.  The probability of unfolding 
%      (or refolding) in the interval <F(i)-dF,F(i+dF)< is pd(i)*dF
%  n: nUmber of rows in force
%  F: Force vector (useful if F is not specified in the call)

  if isempty(force)  % Empty cluster
    pd = NaN;
    F = NaN;
    n = 0;
    return
  end
  if nargin < 3
    plotting = 0;
  end
  edges = floor(min(force)-0.5*dF)+0.5*dF:dF:ceil(max(force) + 0.5*dF)-0.5*dF;
  F = (edges(1:end-1)+edges(2:end))'/2;
  Values = histcounts(force,edges);
  n = sum(Values);
  pd = Values'/n/dF;  % Probability density
  if plotting
    figure;
    histogram(force,edges);
    cla;
    b = bar(F,pd,1); 
    set(b,'facecolor',[0.2,.7,.9])
    title(sprintf('%d events',n))
    ylabel('Probability density')
    xlabel('Force (pN)')
  end
end
