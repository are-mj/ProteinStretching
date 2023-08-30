function p = refold_probability_model(theta,F,T,Fdotmean)
% Refold probaility density using the Bell model
% Input:
%  theta : [dx;log10(k0)];
%  F     : Force vector (pN)
%  T     : Temperature (°C)
%  Fdotmean  : Observed rate of change for force

  dx = theta(1);
  log10k0 = theta(2);

  kB = 0.01380649; % Boltzmann constant, unit: 1e-21 J = 1zJ/K
  beta = 1./(kB*(T+273.15));
  k0 = 10^log10k0;

  % refolding rate:
  kR = k0*exp(-beta*F*dx);
  p = -kR/Fdotmean.*exp(kR/(beta*Fdotmean*dx));
end
