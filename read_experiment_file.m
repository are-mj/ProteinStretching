function [t,f,xx,T,Filename] = read_experiment_file(file)
% Reads in a text file from molecular tweezers protein streching experient
% Input: file
%    Normally on the form <date>/xX.txt'. The full path will be supplied
%      by datafolder.m.
%    May also include the full path.  This will override datafolder.m
% Output:
%    t   : Time (s)  Calculated as CycleCounts*4000;
%    f   : Force (pN)
%    xx  : Two column extent array from A_dist_y and B_dist_y (nm)
%    T   : Temperature if this is coded in the Status column
%    Filename:  Short form <date>'/xX.txt'. 

  % Version 2023-05-07: Allows more flexible file names.  Use full path if
  %    filenames do not conform to <datafolder>/<date>/<xx>.txt
  % Version 2023-01-23: Translate all backslashes in file path to slashes
  %   for compatibility with Unix amd Mac.
  % Version 2: 2023-02-13: reads Columnns A_dist_Y and B_dist_y 
  %   into the two-column array xx.
  % Version 3: 2023-08-15
  
  cps = 4000;  % CycleCounts per second

  file = strrep(file,'\','/');
  if numel(regexp(file,'\/'))<2  % Short form of file 
    file_full = fullfile(datafolder,file);
    file_full = strrep(char(file_full),'\','/');
  else
    file_full = file;
  end
  
  filename_slashes = regexp(file_full,'\/');  % Position of '/' in string file
  Filename = string(file_full(filename_slashes(end-1)+1:end));  % short file name

  % Read file
  fid = fopen(file_full);
  if fid == -1
      error('File %s not found',file)
  end

  textscan(fid,'%s',1,Delimiter='\n');  %skip first line
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
  xx = -[data{xA_column},data{xB_column}];  % Read both trap position columns
  t =  counts/cps;   % Time from experiment start (seconds)
  status = data{statuscolumn};

  T_index = floor(rem(status,1000)/10);  % Number from digits 2 and 3
  T_index(status<1000) = NaN;
  T = NaN(size(T_index));

  % The coding scheme for temperature is given by temperature_code:
  Tlist = temperature_code(Filename);
  if isempty(Tlist)
    T = T_from_COM(file)*ones(size(f));
  else
    for ii = 1:size(Tlist,1)
      T(T_index==Tlist(ii,1)) = Tlist(ii,2);
    end
  end
end

