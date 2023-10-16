function[t,f,xx,T,Filename] = read_experiment_file(file) 
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
  % Version 3: 2023-10-15: 
  %     Automatic detection of large frequency changes
  %     Headers may be found in either line 1 or line 2
  %     Default temp. 20°C if neither Tlist nor COM file available. 
  %     Reads time column (header: 'time(sec)') if it exists.

  
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

  % The headers are normally found in line 1 or 2 of the file:
  headercell = textscan(fid,'%s',1,Delimiter='\n');
  headerline = headercell{1}{1};
  if ~contains(headerline,'X_force')
    % Line 1 does not contanin headers, Try line 2:
    headercell = textscan(fid,'%s',1,Delimiter='\n');
    headerline = headercell{1}{1};   
  end
  if ~contains(headerline,'X_force')
    error('Could not find headers in lines 1 or 2 of %s',file_full)
  end
  headers = regexp(headerline,'\t','split');

  frewind(fid)  
  % read data from line 4 (col 1 in line 3 is often bad)
  data = textscan(fid,'%f %f %f %f %f %f %f %f%*[^\n]','headerlines',3);
  fclose(fid);

  % Define data (column) vectors
  timecol = contains(headers,'time(sec)');
  if any(timecol)
    t = data{timecol};
  else
    countscol = contains(headers,'CycleCount');
    t = data{countscol}*cps;
  end
  forcecol = contains(headers,'Y_force');
  f = -data{forcecol};
  xA_column = contains(headers,'A_dist-Y');
  xB_column = contains(headers,'B_dist-Y');
  xx = -[data{xA_column},data{xB_column}];  % Read both trap position columns
  statuscolumn = contains(headers,'Status');
  status = data{statuscolumn};

  % Try reading temperature 
  Tlist = temperature_code(Filename);
  if ~isempty(Tlist)
    % Read temperature index from digits 2 and 3 in the status column
    T_index = floor(rem(status,1000)/10);  % Number from digits 2 and 3
    T_index(status<1000) = NaN;
    T = NaN(size(T_index));    
    for ii = 1:size(Tlist,1)
      T(T_index==Tlist(ii,1)) = Tlist(ii,2);
    end
  else   
    try
      T = T_from_COM(file)*ones(size(f));
    catch
      warning('Unable to read temperatures. Uses default 20°C');
      T = 20*ones(size(f));
    end
end
