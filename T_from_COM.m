function T = T_from_COM(file)
% Read temperature from the COM file corresponding to the experiment file
%   Backup solution if the temperature is not coded in the staus variable
%   2023-10-16: Modified the construction of COM file name.  Now works for
%   abA.txt (but no longer for aAB.txt)
  file = char(file);
  slashes = regexp(file,'[\/\\]');
  if numel(slashes)<2  % Short version of filename
    file_full = fullfile(datafolder,file);  % datafolder.m must return the folder
    slashes = regexp(file_full,'[\/\\]');
  else
    file_full = file;
  end
  strrep(file_full,'\','/');
  
  fiber = file_full(slashes(end)+1:end-5);
  COMfile = [file_full(1:slashes(end)),fiber,'COM.txt'];
  
  fid = fopen(COMfile);
  c = textscan(fid,'%s','delimiter','\n');
  fclose(fid);
  lines = c{1};
  nlines = numel(lines);
  TA = [];
  for j = 1:nlines
    if contains(lines{j},'temperatureA')
      equalpos = findstr(lines{j},'=');
      numbers = sscanf(lines{j}(equalpos(1)+1:end),'%f');
      TA = [TA;numbers(1)];
      % data = sscanf(lines{j},'%d %*s %*s %f');
      % TA = [TA;data(2)];
    end
  end
  T = round(mean(TA,'omitnan'));
  if isnan(T)
    T = 20;   % Just to choose a default
  end
end