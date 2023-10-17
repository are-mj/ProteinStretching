function [resnorm,theta,theta_std] = fitBell(F,pd,n,Tmean,theta0,Fdotmean)
% Fitting Bell_probability_model to observed probability densities
% Input:
%   F        column array of force values
%   pd       column array of experiment probability densities
%   T        Temperature (°C)
%   theta0   Initial guess for model parameters [dx;log10(k0)]
%   par      parameter struct with fields Fdotmean, nu, model
% F, pd and Fdotmean may be caclulated by probability_density.m
% Output:
%   resnorm  norm of difference between input and model pd
%   theta    Fitted parameters
%   theta_std  Standard deviations for the parameter estimate 
%
% Similar to fit_refold parameters, but returns resnorm as first output
% Useful when optimizing clustering

  if n < 4            % Trick to remove all events fom nearly empty cluster
    resnorm = 10*n;   % May not be robust in all cases
    theta = NaN;
    theta_std = NaN;
    return
  end
  % Optimization parameters:
  opt = optimoptions('lsqcurvefit');
  opt.MaxFunctionEvaluations = 2000;
  opt.Display = 'off';

  lb = [0;-15];
  ub = [500;10000];

  probfun = @(theta,F) Bell_probability_model(theta,F,Tmean,Fdotmean);
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
