function [Tu,Tr,f0,x0,time0,tpeaks,tvalleys] = analyse_file(file,plotting,printing)
% Identify unfolding and refolding events from measurement text file
%   Replaces identify_evnts.m
% Inputs:
%  file: short form (Format '20220718/cA.txt' or '04082022/uA.txt')
%            The path is provided by datafolder.m.
%            Alternatively, use the full file path.   
%  plotting: 1: show plot of force vs. counts, with transition points
%              marked in red
%            0: No plotting (default) 
%  printing: 1: Print summary per file and part (default)
%            0: No printing
% Output:
%  Tu,Tr: Matlab tables of Time, Deltax, Force etc. for unfold and refold
%  f0, x0, time0:  time series of force, trap position and time
%  tpeaks:       Times for detected force maxima. 
%                  Useful for finding refold/unfold pairs
%  tvalleys:     Times for detected force minima. 
%                  Useful for finding unfold/refold pairs
%
% Version 1.1 2023-02-25
% Are Mjaavatten (mjaavatt@gmail.com)
% Version 1.2 2023-03-24: Added error message if called with no arguments
% Version 1.3 2023-07-06: Added column 'Pullingspeed' to Tu, Tr
% Version 1.4 2023-08_26: Removed upper limits for peak width and distance
%                         Added columns Trapx, dFdx and dt to Tu and Tr. 
% Version 1.5 2023-09-04: Automated finding large signal frequency changes
% Version 1.6 2023-09-18: Added output tpeaks. Reduced MinPeakWidth by 17%
% Version 1.7 2023-11-28: Added output tvalleys.f0

  if nargin < 1
    error('Missing input argument: file')
  end
  if nargin <2
    plotting = 0;  % No plotting unless specified
  end    
  if nargin < 3
    printing = 1;
  end

  % Note: 'file' may contain full path or only the short form
  %       'Filename' is always the short form
  [time0,f0,xx0,T,Filename] = read_experiment_file(file);  
  x0 = mean(xx0,2);

  % Define index ranges with similar x oscillation frequency
  % First check for predefined ranges
  r = valid_data_ranges(Filename);
  if r == 0  % No predefined ranges
    % Automatically find changes in extent oscillation frequency
    mx = median(x0)-x0;  % Center
    c = find(mx(2:end).*mx(1:end-1) < 0);  % zero crossing indices
    c(diff(c)<50) = [];  % Remove short crossing intervals caused by noise
    % Find large changes in crossing intervals:
    r = [1;c(ischange(log10(diff(c)),'threshold',1));length(x0)];  
  end
  n_parts = numel(r)-1;

  Tu = cell2table(cell(0,12),'VariableNames',{'Filename','Time','Deltax','Force','Forceshift','Trapx','Fdot','Pullingspeed','Temperature','dFdx','dt','Lineno'});
  Tr = cell2table(cell(0,12),'VariableNames',{'Filename','Time','Deltax','Force','Forceshift','Trapx','Fdot','Pullingspeed','Temperature','dFdx','dt','Lineno'});
  N_traces = 0;
  tpeaks = [];
  tvalleys = [];

  if plotting 
    % Show full time series in grey color as background:
    figure;  
    plot(time0,f0,"Color",.7*[1 1 1]); 
  end

  for part = 1:n_parts
    % Define variables for this part
    % index_range = r(part,1):r(part,2);
    index_range = r(part):r(part+1);
    f = f0(index_range);  % Force
    if max(f)-min(f) < 8
      continue
    end
    x = x0(index_range);  % Extent
    time = time0(index_range);

    % find the dominant frequency of the f series to specify parameters for
    % the findpeaks function
    nu = dominant_frequency((0:length(f)-1),(f-mean(f)));
    if isempty(nu)
      continue
    end
    peak_distance = 1/nu;
    % MinPeakDistance = min(500,fix(peak_distance/2));
    % MinPeakWidth = min(200,fix(peak_distance/5));
    % Fails for slow pulling speeds!
    MinPeakDistance = fix(peak_distance/2);
    MinPeakWidth = fix(peak_distance/6);    
  
    % Force peaks and valleys separatechange_force unfolding and refolding streches:
    [hival,hipos] = findpeaks(f,'MinPeakDistance',MinPeakDistance,'MinPeakWidth',MinPeakWidth);
    [loval,lopos] = findpeaks(-f,'MinPeakDistance',MinPeakDistance,'MinPeakWidth',MinPeakWidth);
    loval = -loval;
    tpeaks = [tpeaks;time(hipos)];
    tvalleys = [tvalleys;time(lopos)];

    pos = sort([hipos;lopos]); % Indices separating stretching and relaxing traces 
    N = numel(pos)-1;     % Number of separate traces
    N_traces = N_traces + N;
    N_unfold = 0;N_refold = 0;
    for i = 1:N
      trace = pos(i):pos(i+1);
      if numel(trace)<100
          continue;
      end
      s.f = f(trace);
      % The recorded x is negative. Transform so that x increases with f:
      s.x = max(x(trace))-x(trace);
      s.t = time(trace);
      dt = s.t(end)-s.t(1);
      % Discard empty or very short or very long traces:
      if isempty(s.f) || dt < 0.1 || dt > 35  
        continue
      end
      [k1,Force,Deltax,Fdot,Forceshift,Pullingspeed,dFdx,dt] = analyse_trace(s); 
      if k1 < 1
        continue;  % No transition found
      end
      Trapx = s.x(k1);
      index = k1 + trace(1)-1;  
      Time = time(index);
      Lineno = index_range(index)+4;   % Line number in file
      if numel(T) > 1
        Temperature = T(Lineno);
      else
        Temperature = T;
      end
      if sign(s.f(end)-s.f(1)) == 1
        Tu = [Tu;table(Filename,Time,Deltax,Force,Forceshift,Trapx,...
          Fdot,Pullingspeed,Temperature,dFdx,dt,Lineno)];
        N_unfold = N_unfold + 1;
      else
        Tr = [Tr;table(Filename,Time,Deltax,Force,Forceshift,Trapx, ...
          Fdot,Pullingspeed,Temperature,dFdx,dt,Lineno)];  
        N_refold = N_refold + 1;
      end
    end
    ratio = (N_unfold+N_refold)/(numel(pos)-1);
    if printing
      fprintf('Part %4d: traces: %4d  Unfoldings %4d  Refoldings %4d, Ratio: %4.2f\n',...
        part,numel(pos)-1,N_unfold,N_refold,ratio)
    end

    if plotting
      hold on;
      plot(time,f,'b')
      h = plot(Tu.Time,Tu.Force,'or');
      set(h,'MarkerFaceColor','r','MarkerSize',4);
      plot(Tr.Time,Tr.Force,'ok');     
      h = plot(time([hipos;lopos]),[hival;loval],'ok');
      set(h,'MarkerFaceColor','k','MarkerSize',2);   
      hold off
    end 
  end  % end part
  
  N_unfold = size(Tu,1);
  N_refold = size(Tr,1);
  ratio = (N_unfold+N_refold)/N_traces;
  if printing
    fprintf('%s: Found %4d unfoldings and %4d refoldings in %4d traces. Ratio: %5.2f\n', ...
      Filename,N_unfold,N_refold,N_traces,ratio);
  end
  if plotting 
    ylim([0,55]);  % Common force range for all plots
    titlestr = sprintf(Filename);
    title(titlestr,'interpreter','none') 
    ylabel('Force (pN)');
    xlabel('Time (s)')
  end  
end

