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