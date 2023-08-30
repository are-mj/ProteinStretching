function DQcal = convert_energy(DQ)
% Convert energy unit 
%   from: zJ/molecule (zJ = 1e-21J = nm*pNJ)
%   to:   kcal/mol
  % NA = 6.02214e23; % molecules/mol
  % kcal = 4184; % J
  % nm = 1e-9;   % m
  % pN = 1e-12;  % N
  % zJ = nm*pN;  % J
  % conversion = NA*zJ/kcal; % zJ to kcal/mol;
  conversion = 0.1439326;
  DQcal = conversion*DQ;
end