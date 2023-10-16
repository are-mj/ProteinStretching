function x = wlc_inverse(f,P,T,L0)
% Worm-like chain extension as function of the force
% Input:
%  F: Force (pN)
%  P: Persistence length (nm)
%  T: Temperature (K)
%  L0: Contour length (Molecule length at maximun extension) (nm)
% Output: 
%  x: Extension of molecule (nm) 

  x = zeros(size(f));
  for i = 1:length(f)
    fun = @(x) wlc(x,P,T,L0)-f(i);
    x(i) = fzero(fun,10);
  end
end