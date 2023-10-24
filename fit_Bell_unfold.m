function [theta,theta_std,resnorm] = fit_Bell_unfold(F,pd,Tmean,Fdot,theta0)
% Fitting Bell_unfold_probability to observed probability densities
% Input:
%   F        column array of force values
%   pd       column array of experiment probability densities
%   Tmean    Temperature (°C)
%   Fdot     mean value of dF/dt before unfolding 
%   theta0   Initial guess for model parameters [dx;log10(k0)]
% F, pd, Tmean, and Fdot may be caclulated by probability_density.m
% Output:
%   theta     Fitted parameters
%   theta_std Standard deviations for the parameter estimate 
%   resnorm   Norm of difference between input and model pd

  % Optimization parameters:
  opt = optimoptions('lsqcurvefit');
  opt.MaxFunctionEvaluations = 2000;
  opt.Display = 'off';

  lb = [0;-15];
  ub = [500;10000];

  probfun = @(theta,F) Bell_unfold_probability(theta,F,Tmean,Fdot);
  [theta,resnorm,resid,exitflag,~,~,J] = ...
    lsqcurvefit(probfun,theta0,F,pd,lb,ub,opt);
  if exitflag < 1
    error('lsqcurvefit problems. Exitflag: %d',exitflag)
  end 
  if nargout > 2
    ci = nlparci(theta,resid,'jacobian',J);  % 95% confidence interval
    theta_std = (ci(:,2)-theta)/fzero(@(x)normcdf(x)-0.975,-1); % Standard deviation
  end
end
