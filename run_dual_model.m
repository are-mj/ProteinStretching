% run_dual_model.m

load Top7_Tables.mat Tun Tre
speed = Tun.Pullingspeed;
Temperature = Tun.Temperature;
fast = speed>400;
med = speed<400 & speed>40;
slow = speed<40;
clear cases;
cases(7) = struct();  % Initiaiise
cases(1).selected = fast&Temperature<8 & Temperature>=4;
cases(1).speedtext = "Fast";
cases(1).temptext = "Low";
cases(2).selected = fast&Temperature>18 & Temperature < 22;  
cases(2).speedtext = "Fast";
cases(2).temptext = "High";
cases(3).selected = med&Temperature<8;
cases(3).speedtext = "Normal";
cases(3).temptext = "Low";
cases(4).selected = med&Temperature>8 & Temperature < 15;
cases(4).speedtext = "Normal";
cases(4).temptext = "Medium";
cases(5).selected = med&Temperature>18 & Temperature < 22;
cases(5).speedtext = "Normal";
cases(5).temptext = "High";
cases(6).selected = slow&Temperature<8;
cases(6).speedtext = "Slow";
cases(6).temptext = "Low";
cases(7).selected = slow&Temperature>18; 
cases(7).speedtext = "Slow";
cases(7).temptext = "High";

% Output table heading:
fprintf('  Speed    Temp. Cluster Model     ΔG‡           x‡          log10(k0)       n    Residual norm\n');
splitforce = zeros(7,1);
outbell = zeros(14,6);
outdudko = zeros(14,8);
for i = 7
  % fit_dual(Tun,cases(i));
  fit_dual(Tre,cases(i));
end
