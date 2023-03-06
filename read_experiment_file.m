function [file,t,f,xx,T] = read_experiment_file(filename)
% Reads in a text file from molecular tweezers protein streching experient
% Input: filename
%    Normally on the form '04082022/uA.txt'. The full path will be supplied
%      by datafolder.m.
%    May also include the full path.  This will override datafolder.m
% Output:
%    file:  Format 'mmddyyyy/<xx>.txt'
%       Alternatively, from the data folder part of filename.  
%    t   : Time (s)  Calculated as CycleCounts*4000;
%    f   : Force (pN)
%    xx  : Two column extent array from A_dist_y and B_dist_y (nm)
%    T   : Temperauser if this is coded in the Status column

  % Version 2023-01-23: Translate all backslashes in file path to slashes
  %   for compatibility with UNix amd Mac.
  % Version 2: 2023-02-13: reads Columnns A_dist_Y and B_dist_y 
  %   into the two-column array xx.
  
  cps = 4000;  % CycleCounts per second

  file = strrep(filename,'\','/');
  if numel(regexp(file,'\/'))<2  % Short form of filename ('<date>/xx.txt')
    file_full = strrep(char(fullfile(datafolder,filename)),'\','/');
    filename_slashes = regexp(file_full,'\/');  % Position of '/' in string file
    file = string(file_full(filename_slashes(end-1)+1:end));
  end

  % Read file
  fid = fopen(file);
  if fid == -1
      error('File %s not found',file)
  end
  % First line gives experiment date. THis is no longer used
  try
    dateline = textscan(fid,'%s',1,Delimiter='\n');  
    slashes = cell2mat(regexp(dateline{1},'/'));
    datestr = dateline{1}{1}(slashes(1)-2:slashes(2)+2);
    date = datetime(datestr,'inputformat','MM/dd/yy');
  catch
    try 
     parts = strsplit(file,'/');
     date = datetime(parts{1},'inputformat','MMddyyyy');
    catch
      date = datetaime('01-Jan-2020');  % Unknown date
    end
  end
  % Line 2 gives column headers
  headercell = textscan(fid,'%s',1,Delimiter='\n');
  headerline = headercell{1}{1};
  headers = regexp(headerline,'\t','split');
  xA_column = contains(headers,'A_dist-Y');
  xB_column = contains(headers,'B_dist-Y');
  statuscolumn = contains(headers,'Status');
  frewind(fid)

  % read data from line 4 (col 1 in line 3 is often bad)
  data = textscan(fid,'%f %f %f %f %f %f %f %f%*[^\n]','headerlines',3);
  fclose(fid);

  % Define data (column) vectors
  counts = data{1};
  f = -data{3};      % Force
  xx = -[data{xA_column},data{xB_column}];  % Read both extent columns
  t =  counts/cps;   % Time from experiment start (seconds)
  status = data{statuscolumn};

  T_index = floor(rem(status,1000)/10);  % Number from digits 2 and 3
  
% The coding scheme for temperature is not always the same
% More Tlist definitions may be needed in the future
  if date == datetime('18-Jul-2022') 
    Tlist = [0,4.4;2,9.88;3,12.12;4,14.89;6,18.89;16,31.7];
  else
    Tlist = [(0:2:12)',[4;6;10;14;16;20;23]];
  end
  T = nan(size(T_index));
  for ii = 1:size(Tlist,1)
    T(T_index==Tlist(ii,1)) = Tlist(ii,2);
  end
end