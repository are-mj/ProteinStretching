function s = extract_trace(filename,t,plotting)
% Extract one tretch containing time t from measurement text file
% Inputs:
%  file:     data file <xx>.txt.  Include the full path if the file is not
%            in the current folder
%  t     Time for a point within the desired trace
% Output:
%  s: struct containing data for the selected trace
%    s.r:  Index range
%    s.f:  Force
%    s.x:  Distance
%    s.t:  time
%    s.file: filename
%

% Version 0.1, 2022-09-13
% Version 0.2, 2022-09-20
% Version 1.0, 2022-10-06
% Version 1.1, 2023-08-15  % removed lower limit for peak spacing
% Are Mjaavatten 

  file = strrep(filename,'/','\');
  if numel(regexp(file,'\\'))<2  % Short form of filename ('<date>/xx.txt')
    file = fullfile(datafolder,file);  % Full file name
  end

  if nargin < 3
    plotting = 0;
  end

  file = strrep(filename,'\','/');
  if numel(regexp(file,'\/'))<2  % Short form of filename ('<date>/xx.txt')
    file = char(fullfile(datafolder,filename));
    file = strrep(file,'\','/');
  end
  [t0,f0,xx0,T] = read_experiment_file(file);
  x0 = mean(xx0,2);
  if t < t0(1) || t > t0(end)
    error('Requested time (%.2f) is outside the file range %.2f - %.2f',t,t0(1),t0(end));
  end
  if t0(1)>t0(2)
    f0(1) = [];
    x0(1) = [];
    t0(1) = [];
  end
  ny0 = numel(f0);  


  % Read parameters for file
  r = valid_data_ranges(file);  
  if isempty(r)
    r = [1,numel(f0)];
  end  

  index0 = interp1(t0,1:ny0,t,'nearest');
  part = find(r(:,1)<index0,1,'last');
  index = index0 - r(part,1) +1;  % Local index within part

  f = f0(r(part,1):r(part,2));  % Force
  x = x0(r(part,1):r(part,2));  % Extent
  lines = r(part,1):r(part,2);   % Line numbers in file
  time = t0(lines);
  nu = dominant_frequency((0:length(f)-1),(f-mean(f)));
  peak_distance = 1/nu;
  % MinPeakDistance = min(500,fix(peak_distance/2));
  % MinPeakWidth = min(200,fix(peak_distance/5));
  % Fails for slow pulling speeds!
  MinPeakDistance = fix(peak_distance/2);
  MinPeakWidth = fix(peak_distance/5);  

  % Locate force peaks and valleys
  [~,hipos] = findpeaks(f,'MinPeakDistance',MinPeakDistance,'MinPeakWidth',MinPeakWidth);
  [~,lopos] = findpeaks(-f,'MinPeakDistance',MinPeakDistance,'MinPeakWidth',MinPeakWidth);
  pos = sort([hipos;lopos]); % Indices of switch point between fold type 
  N = find(index < pos,1);
  if isempty(N) | N < 2
    s = NaN;
    return
  end
  trace = pos(N-1):pos(N);
  if numel(trace)<10
    error('No menaingful trace found')
  end
  s.f = f(trace);
  s.x = max(x(trace))-x(trace);
  s.t = time(trace);
  s.r = trace;
  k = analyse_trace(s);
  s.k = k;
  s.file = filename;
  if plotting
    figure;
    subplot(211);plot(s.t,s.f);
    xlabel('time (s)')
    ylabel('Force (pN)')
    subplot(212);plot(s.f)
    xlabel('index')
    ylabel('Force (pN)')
    if k>0
      subplot(211),hold on;
      plot(s.t(k),s.f(k),'.r','MarkerSize',10);
      subplot(212),hold on;
      plot(k,s.f(k),'.r','MarkerSize',10); 
    end
  end
end
