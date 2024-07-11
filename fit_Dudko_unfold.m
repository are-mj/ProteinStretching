function [theta,theta_std,resnorm,resid] = fit_Dudko_unfold(F,pd,T,Fdot,theta0,par)
% Fitting parameters for the DHS or CHS unfolding model
%   DHS: Dudko, Hummer & Szabo (2006), CHS: Cossio, Hummer & Szabo (2016)
% Input:
%   F        array of force values
%   pd       array of experiment probability densities
%   T        Temperature (°C)
%   Fdot     Mean value of df/dt before unfolding
%   theta0   Initial guess for model parameters [DG;dx;log10(k0)]
%   par      parameter struct with fields Fdot, nu, model
%     nu:  1/2: quadratic potentiial, 2/3 cubic potential
%     model alternatives:
%       'DHS' (Dodko, Hummer & Szabo, 2006)  
%       'CHS' (Cossio, Hummer & Szabo, 2016) 
%   F, pd, T and Fdot may be calculated by probability_density.m
% Output:
%   theta     Fitted parameters
%   theta_std Standard deviations for the parameter estimates 
%   resnorm   Norm of difference between input and model pd

% Version 4.1,October 2023

  % Optimization parameters:
  opt = optimoptions('lsqcurvefit');
  opt.MaxFunctionEvaluations = 2000;
  opt.Display = 'off';
  
  % theta0: initial guess for theta = [DG;dx;log10k0]
  DG_0 = theta0(1);
  dx_0 = theta0(2);
  log10k0_0 = theta0(3);

  a_0 = par.nu*dx_0/DG_0;  % Use a as parameter.  Setting limits on a
                           % lets us avoid complex numbers in kfun.
  % Alternative parameter vector for use in fminsearch:
  thetacalc0 = [DG_0;a_0;log10k0_0];
  % Symmetric error bars in log10(k) avoids nonsensical error bars
  
  % Lower and upper parameter bounds
  lb = [0;0;-25];  %If theta(3) = log10(k)
  ub = [500;0.9999/F(end);1];  % Make sure that 1-a*F > 0 in the probalility expression
  
  % Function for calculating model probabilities
  probfun = @(thetacalc,F)Dudko_unfold_probability(thetacalc,F,T,Fdot,par);
  
  % Fit model parameters to data:
  [thetacalc,resnorm,resid,exitflag,~,~,J] = ...
    lsqcurvefit(probfun,thetacalc0,F,pd,lb,ub,opt);
  if exitflag < 1
    % warning('lsqcurvefit problems. Exitflag: %d',exitflag)
    theta = NaN;
    theta_std = NaN;
    resnorm = NaN;
    Fplot = NaN;
    pdplot = NaN;
    return
  end
  DG = thetacalc(1);
  a = thetacalc(2);
  log10k0 = thetacalc(3);
  dx  = a*DG/par.nu;
  theta = [DG;dx;log10k0];
  
  % ***  Confidence interval for identified perameters ***
  % J is the Jacbian of f(thetacalc) = taufun(thetacalc,...).  
  % JG is the Jacbian of fG(theta) = f(thetacalc(theta))
  da_dDG = -par.nu*dx/DG^2;
  da_dx = par.nu/DG;
  JG = [J(:,1)+J(:,2)*da_dDG,J(:,2)*da_dx,J(:,3)];
  ci = nlparci(theta,resid,'jacobian',JG);  % 95% confidence interval
  theta_std = (ci(:,2)-theta)/fzero(@(x)normcdf(x)-0.975,-1); % Standard deviation
  Fplot = linspace(0,round(F(end)/5)*5+5);
  pdplot = probfun(thetacalc,Fplot);  % Model output with fitted paramters
end