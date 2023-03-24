function [Tu,Tr,f0,x0,time0] = analyse_file(filename,plotting)
% Identify unfolding and refolding events from measurement text file
%   Replaces identify_evnts.m
% Inputs:
%  file:     data file <xx>.txt.  Include the full path if the file is not
%            in the current folder
%  plotting: 1: show plot of force vs. counts, with transition points
%              marked in ref
%            0: No plotting (default) 
% Output:
%  Tu,Tr: Matlab tables of Time, Deltax, Force etc. for unfold and refold
%  f0, xo, time0:  time series of force extent and time
%
% Version 1.1 2023-02-25
% Are Mjaavatten (are@mjaavatten.com)
% Version 1.2 2023-03-24: Added error message if called with no arguments

  if nargin < 1
    error('Missing input argument: filename')
  end
  if nargin <2
    plotting = 0;  % No plotting unless specified
  end    

  [Filename,time0,f0,xx0,T] = read_experiment_file(filename);  
  x0 = mean(xx0,2);
  
  % Read valid index ranges for file if they are defined
  r = valid_data_ranges(Filename);
  if isempty(r)
    r = [1,numel(f0)];
  end 
  n_parts = size(r,1);  % Number of parts to be analysed

  Tu = cell2table(cell(0,8),'VariableNames',{'Filename','Time','Deltax','Force','Forceshift','Fdot','Temperature','Lineno'});
  Tr = cell2table(cell(0,8),'VariableNames',{'Filename','Time','Deltax','Force','Forceshift','Fdot','Temperature','Lineno'});

  if plotting 
    % Show full time series in grey color as background:
    figure;  
    plot(time0,f0,"Color",.7*[1 1 1]); 
  end

  for part = 1:n_parts
    % Define variables for this part
    index_range = r(part,1):r(part,2);
    f = f0(index_range);  % Force
    x = x0(index_range);  % Extent
    time = time0(index_range);

    % find the dominant frequency of the f series to specify paramters for
    % the findpeaks function
    nu = dominant_frequency((0:length(f)-1),(f-mean(f)));
    peak_distance = 1/nu;
    MinPeakDistance = min(500,fix(peak_distance/2));
    MinPeakWidth = min(200,fix(peak_distance/5));
  
    % Force peaks and valleys separatechange_force unfolding and refolding streches:
    [hival,hipos] = findpeaks(f,'MinPeakDistance',MinPeakDistance,'MinPeakWidth',MinPeakWidth);
    [loval,lopos] = findpeaks(-f,'MinPeakDistance',MinPeakDistance,'MinPeakWidth',MinPeakWidth);
    loval = -loval;

    pos = sort([hipos;lopos]); % Indices separating stretching and relaxing stretches 
    N = numel(pos)-1;     % Number of separate stretches
    
    N_unfold = 0;N_refold = 0;
    for i = 1:N
      stretch = pos(i):pos(i+1);
      n_stretch = numel(stretch);
      if n_stretch >10000 || n_stretch < 10  % Unrealstic stretch
        continue;  
      end
      s.f = f(stretch);
      % The recorded x is negative. Transform so that x increases with f:
      s.x = max(x(stretch))-x(stretch);
      s.t = time(stretch);
      if isempty(s.f)
        continue
      end
      [k1,Force,Deltax,Fdot,Forceshift] = analyse_stretch(s); 
      if k1 < 1
        continue;
      end
      index = k1 + stretch(1)-1;  
      Time = time(index);
      Lineno = index_range(index)+4;   % Line number in file
      Temperature = T(Lineno);
      if sign(s.f(end)-s.f(1)) == 1
        Tu = [Tu;table(Filename,Time,Deltax,Force,Forceshift,Fdot,Temperature,Lineno)];
        N_unfold = N_unfold + 1;
      else
        Tr = [Tr;table(Filename,Time,Deltax,Force,Forceshift,Fdot,Temperature,Lineno)];  
        N_refold = N_refold + 1;
      end
    end
%     ranges(part,2) = size(Tu,1)+1;
%     ranges(part,4) = size(Tr,1)+1;
    ratio = (N_unfold+N_refold)/(numel(pos)-1);
    fprintf('Part %d: Stretches: %d  Unfoldings %d  Refoldings %d, Ratio: %4.2f\n',...
      part,numel(pos)-1,N_unfold,N_refold,ratio)

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
%     allhi = [allhi;hipos];
  end
  
  N_unfold = size(Tu,1);
  N_refold = size(Tr,1);
  ratio = (N_unfold+N_refold)/(numel(pos)-1);
  fprintf('%s: Found %4d unfoldings and %4d refoldings in %4d stretches. Ratio: %5.2f\n', ...
    Filename,N_unfold,N_refold,numel(pos)-1,ratio);
  if plotting 
    titlestr = sprintf(Filename);
    title(titlestr,'interpreter','none') 
    ylabel('Force (pN)');
    xlabel('Time (s)')
  end  
end

