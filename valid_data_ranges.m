function r = valid_data_ranges(file)
% Start and end indices for the parts of input data to be analysed.
% Input: short file name <date>/xX.txt.  NOTE: date format may vary
% Output: r (n by 2 array of start end end indices for parts)
%
%   data index normally corresponds to line number in the file minus 4.
%   r: n by 2 array
%
% Used by analyse_file.m
%
% To specify parts in a new file, use read_experiment_file.m 
%   and plot f.  Use the zoom and data tip tools to specify start
%   and end indices of parts with relatively homogeneous force osclillation
%   frequency. If no r is specified for a given file, the whole 
%   data set is used by analyse_file

  file = char(file);  % Convert if string
  file = strrep(file,'\','/');
  s = strfind(file,'/');
  if numel(s)>1
    file = file(s(end-1)+1:end);  % Short form of file name
  end

  r = 0;

  switch file
    case '02142022/bA.txt'
      r = [1;425500;470850;531810];
    case '07012022/zA.txt'
      r = [1;75350;111910;390500;455600;528000;574200;807103];
    case '07142022/aA.txt'
      r = [1;197200;269025;355000];
  end
end