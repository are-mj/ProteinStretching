% run_all.m
% analyse all experiment results
plotting = 0;
files = allfiles;

Tun = [];Tre = [];
for i = 1:numel(files)
  [Tu,Tr,N] = analyse_file(files{i},plotting);
  Tun = [Tun;Tu];
  Tre = [Tre;Tr];
end
