function [Tun,Tre,pairs] = analyse_many(files,plotting)
% Analyse a set of experiment data files into tables Tun and Tre
% pairs: Row numbers in Tre and Tun for event pairs from connected
%   relax-stretch traces
%
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

function pairs = findpairs(tpeaks,Tu,Tr)
% Find refolding and unfolding events in the same relax-stretch cycle
% i.e.: No force peaks between the two events
% All inputs are outputs from analyse_file
% Output pairs is a column array of corresponding indices [i,j,k]
%  Tr(i,:) Renfolding table row
%  Tr(j,:) Unfolding table row
%  tpeaks(k): Time for the following valley
  pairs = [];
  for i = 1:numel(Tr.Time)
    j = find(Tu.Time>Tr.Time(i),1);
    k = find(tpeaks>Tr.Time(i),1);
    next_unfold = Tu.Time(j);
    next_peak = tpeaks(k);
    if next_unfold < next_peak
      pairs = [pairs;i,j,k];
    end
  end
end