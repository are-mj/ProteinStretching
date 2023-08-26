function [Tun,Tre] = analyse_many(files,plotting)
% Analyse a set of experiment data files into tables Tun and Tre
  if nargin < 2
    plotting = 0;
  end
  Tun = [];
  Tre = [];
  for ii = 1:length(files)
    [Tu,Tr] = analyse_file(files{ii},plotting);
    Tun = [Tun;Tu];
    Tre = [Tre;Tr];
    if plotting
      drawnow;
    end
  end
end