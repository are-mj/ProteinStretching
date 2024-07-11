function [t,f,xx,T,shortname] = read_experiment_file(filename) 
% Reads in a text file from molecular tweezers protein streching experient
% Input: filename
%    Normally on the form 'xx.txt' or '<subfolder>/xX.txt'. 
%      The full path will be supplied by datafolder.m.
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
  % Version 3.1: 2023-10-26
  %     textscan format expression can now handle any number of columns

  
  cps = 4000;  % CycleCounts per second

  d = dir(filename);
  if isempty(d)  % File not found in current path
    filename = fullfile(datafolder,filename);
    d = dir(filename);
    if isempty(d)
      error('File %s nor found\n',filename)
    end
  end
  filename = strrep(filename,'\','/');  % Use Unix separator
  filename_slashes = regexp(filename,'\/');  % Position of '/' in files
  fn = char(filename);  % Translate to character array
  if numel(filename_slashes) < 2
    shortname = fn;  % form: 'xx.txt' or '<subfolder>/xX.txt'
  else
    shortname = fn(filename_slashes(end-1)+1:end); % '<subfolder>/xX.txt'
  end
  shortname = string(shortname);  % Convert back to string
																			  																			 
  % Read filename:
  fid = fopen(filename);
  if fid == -1
      error('File %s not found',filename)
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
    error('Could not find headers in lines 1 or 2 of %s',filenamel)
  end
  headers = regexp(headerline,'\t','split');

  frewind(fid)  
  % read data from line 4 (col 1 in line 3 is often bad)
  % data = textscan(fid,'%f %f %f %f %f %f %f %f%*[^\n]','headerlines',3);
  data = textscan(fid,repmat('%f ',1,numel(headers)),'headerlines',3);
  fclose(fid);

  % Define data (column) vectors
  timecol = contains(headers,'time(sec)');
  if any(timecol)
    t = data{timecol};
  else
    countscol = contains(headers,'CycleCount');
    t = data{countscol}/cps;
  end
  forcecol = contains(headers,'Y_force');
  f = -data{forcecol};
  xA_column = contains(headers,'A_dist-Y');
  xB_column = contains(headers,'B_dist-Y');
  xx = -[data{xA_column},data{xB_column}];  % Read both trap position columns
  % Temperature 
  Tbath = T_from_COM(shortname); % Temperature outside cell
  if isempty(Tbath)
    warning('Found no temperature in COM file. Temperatures nor read');
  end
  Tlist = temperature_code(shortname); % Calibrated T when Tbath = Tlist(1,2);
  if isempty(Tlist)
    Tlist = Tbath;
  end
 
  % Read heater setting from digits 2 and 3 in the status column
  statuscolumn = contains(headers,'Status');
  status = data{statuscolumn};   
  heater_setting = floor(rem(status,1000)/10);  % Number from digits 2 and 3
  heater_setting(status<1000) = NaN;

  T = NaN(size(heater_setting)); 
  if numel(Tlist) == 1						   
    T = Tlist*ones(size(heater_setting));
  else
    for ii = 1:size(Tlist,1)
      T(heater_setting==Tlist(ii,1)) = Tbath + Tlist(ii,2) - Tlist(1,2);				   
    end
  end
end
