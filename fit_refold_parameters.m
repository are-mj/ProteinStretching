function [theta,theta_std,resnorm,Fplot,pdplot] = fit_refold_parameters(F,pd,T,theta0,par)
% Fitting parameters for Bell refolding model
% Input:
%   F        column array of force values
%   pd       column array of experiment probability densities
%   T        Temperature (°C)
%   theta0   Initial guess for model parameters [dx;log10(k0)]
%   par      parameter struct with fields Fdotmean, nu, model
% F, pd and Fdotmean may be caclulated by probability_density.m
% Output:
%   theta    Fitted parameters
%   theta_std  Standard deviations for the parameter estimate 
%   resnorm  norm of difference between input and model pd
%   Fplot, pdplot: Data for plotting of fitted model

% Version 1.0 2023-08-13

  % Optimization parameters:
  opt = optimoptions('lsqcurvefit');
  opt.MaxFunctionEvaluations = 2000;
  opt.Display = 'off';

  lb = [0;0];
  ub = [50;10];

  probfun = @(theta,F)refold_probability_model(theta,F,T,par.Fdotmean);

  % Fit model parameters to data:
  [theta,resnorm,resid,exitflag,~,~,J] = ...
    lsqcurvefit(probfun,theta0,F,pd,lb,ub,opt);
  if exitflag < 1
    error('lsqcurvefit problems. Exitflag: %d',exitflag)
  end  
  
    % ***  Confidence interval for identified perameters ***
  ci = nlparci(theta,resid,'jacobian',J);  % 95% confidence interval
  theta_std = (ci(:,2)-theta)/fzero(@(x)normcdf(x)-0.975,-1); % Standard devation

  Fplot = linspace(0,round(F(end)/5)*5+5);
  pdplot = probfun(theta,Fplot);  % Model output with fitted paramters

end
