function p = Bell_unfold_probability(theta,F,T,Fdot)
% Unfold probability density using the Bell model
% Input:
%  theta : [dx;log10(k0)];
%  T     : Temperature (°C)  (scalar)
%  F     : Force vector (pN)
%  Fdot  : Observed rate of change for force

  dx = theta(1);
  log10k0 = theta(2);

  kB = 0.01380649; % Boltzmann constant, unit: 1e-21 J = 1zJ/K
  beta = 1./(kB*(T+273.15));
  k0 = 10^log10k0;

  % Unfolding rate:
  kU = k0*exp(beta*F*dx);
  p = kU/Fdot.*exp(k0/(beta*Fdot*dx)*(1-exp(beta*F*dx)));
end