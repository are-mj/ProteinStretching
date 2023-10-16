function [Tun,Tre,pairs] = analyse_many(files,plotting)
% Analyse a set of experiment data files into tables Tun and Tre
% pairs: Row numbers in Tre and Tun for event pairs from connected
%   release-stretch traces
  if nargin < 2
    plotting = 0;
  end
  Tun = [];
  Tre = [];
  pairs = [];
  for ii = 1:length(files)
    [Tu,Tr,~,~,~,peaks] = analyse_file(files(ii),plotting);
    if nargout > 2
      pp = findpairs(peaks,Tu,Tr);
      if ~isempty(pp)
        pairs = [pairs;pp(:,1:2)+[height(Tre),height(Tun)]];
      end
    end
    Tun = [Tun;Tu];
    Tre = [Tre;Tr];
    if plotting
      drawnow;
    end
  end
end