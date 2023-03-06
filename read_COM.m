function [counts,TA,TB] = read_COM(file)
  if numel(regexp(file,'\/'))<2  % Short version of filename
    file = fullfile(datafolder,file);  % datafolder.m must return the folder
  end
  fid = fopen(file);
  c = textscan(fid,'%s','delimiter','\n');
  fclose(fid);
  lines = c{1};
  nlines = numel(lines);
  counts = [];
  TA = [];
  TB = [];
  for j = 1:nlines
    if contains(lines{j},'temperatureA')
      data = sscanf(lines{j},'%d %*s %*s %f %*s %*s %f');
    if isempty(data)
      nn = regexp(lines{j},'\d?\.?\d+','match');
      data = [NaN,cellfun(@str2num,nn)];
    end
    counts = [counts;data(1)];
    TA = [TA;data(2)];
    TB = [TB;data(3)];
    end
  end
end