function [k,force,dx,Fdot,shift,pullingspeed,dFdx,dt,noise] = analyse_trace(s,plotting,params)
% Find transition point and distance in a unfolding or refolding 
% Inputs;:
%  s:  input struct with fields:
%    x: extent (trap position)
%    f: force 
%    t: time
%  plotting: plot f(x) with linear fits and dx 
%  params: Struct of parameters overriding default values
%     params.multiple = 1: Handle multiple unfoldings and refoldings in
%     same trace
%    
% Outputs
%  k      Index of tansition point in f. k < 1: No transition found
%         Negative k values indicate exit point, for debugging
%  force: Force at unfoldimg or refolding event
%  dx:    distance
%  Fdot:  Force rate of change
%  shift: Force shift at unfolding/refolding transition
%  pullingspeed: mean rate of change of trap position 
%
% Requires the following functions in addition to standard Matlab
%  -findpeaks from the Signal Processing Toolbox
%  -movingslope by John D'Errico:
%    https://www.mathworks.com/matlabcentral/fileexchange/16997-movingslope
%    MATLAB Central File Exchange. Retrieved October 2, 2022.

% Version 2.12
% Author: Are Mjaavatten (mjaavatt@gmail.com)
% Date:   2023-08-30

% Change history
% v 2.8:  Extra output: shift
%         No rounding of recovery_ix.  Avoids spurious discretization of dx
% v 2.9   Restructured to reflect analyse_trace_demo.mlx
% v 2.10  Added ouput: pullingspeed
% v 2.11  Corrected typos.
% v 2.12  Adjusted noise sensitivity

%% Algortithm tuning parameters
  % These parameter values are chosen by trial and error and seem to give
  % reasonable results i most cases
  algpar.minforcerange = 10;  % Discard traces where the force range is small
  algpar.minpointforfit = 3;  % Minimum number of points for linear fit
  algpar.maxsloperatio = 0.5; % Maximum ratio of slopes beore and after transition
  algpar.multiple = 0;        % Modify algorithm to handle multple transitions

  % Transitions are unlikely near the start of a trace, so we exclude
  % potential transitions if abs(f-fstart) < algpar.dfstart
  algpar.dfstart = 4;  
  % At the end of a REFOLD trace noise makes detections difficult:
  algpar.dfend = 1;  % Make sure that abs(f-fend) > algpar.dfend at transition.
  % No end restrictions for unfold traces

  % Low f or x slope parts at start or end of a trace may confuse the  
  % detection algorithm, so eliminate start and end parts where the slope
  % is less that algpar.slopefrac*median(slope):
  algpar.slopefracf = 0.7;
  algpar.slopefracx = 0.8;

%% Handle defaults
  if nargin < 2
    plotting = false;
  elseif isempty(plotting)
    plotting = false;
  end
  if nargin < 3
    params = [];
  end
  if isempty(params)
    params = struct;
  end
  fields = fieldnames(params);
  for ii = 1:numel(fields)
    algpar.(fields{ii}) = params.(fields{ii});
  end

  % Parameters for moving slope calculation and slope peak detection
  % Again, these are not set in stone and may need adjusting
  nf = numel(s.f);
  if nf > 1000
    slopelength = 30;
    minpeak = 0.03;
  elseif nf > 300
    slopelength = 10;
    minpeak = 0.05;
  else
    slopelength = 5;
    minpeak = 0.1;
  end

%% Preprocessing  
  % Make sure all output variables are set:
  force = NaN;
  dx = NaN;
  Fdot = NaN;
  shift = NaN;
  pullingspeed = NaN;
  dFdx = NaN;
  noise = NaN;

% Make sure x>0 and f both increase or decrease together
  s.x = s.x - min(s.x);
  
% Change sign of x and f for refolding traces
% This lets us use the same algorithm for bot unfolding and refolding
  sgn = sign(s.f(end) - s.f(1));
  f = sgn*s.f;
  x = sgn*s.x;
  t = s.t;
  dt = (t(end)-t(1))/(numel(t)-1);

% Do not accept very short traces
  df = f(end)-f(1);  % range of forces
  if(df<algpar.minforcerange)  % Skip if force amplitude is small (probably noise only)
    k = -1;
    return
  end

% Eliminate bad parts at start and end
  slope = movingslope(f,slopelength);
  medslope = median(slope);
  xslope = movingslope(x,slopelength);
  xmedslope = median(xslope);
  startslope = find(slope>algpar.slopefracf*medslope,1);  % End of inital flat part 
  startxslope = find(xslope>algpar.slopefracx*xmedslope,1);
  startmin = find(abs(f(1)-f)>algpar.dfstart,1);
  start = max([1,startslope,startxslope,startmin]);
  stopslope = find(slope>algpar.slopefracf*median(slope),1,'last'); 
  stopxslope = find(xslope>algpar.slopefracx*median(xslope),1,'last'); 
  stopmin = find(f(nf) - f > algpar.dfend,1,'last');
  stop = min([stopmin,stopslope,stopxslope,nf]);
  validrange = start:stop;
  if stop-start < 15
    k = 0;
    return
  end
  pullingspeed = abs((s.x(stop)-s.x(start))/(s.t(stop)-s.t(start)));

%% Find transitions (unfolding or refolding events)

% Find slope changes

  slope = slope(validrange);
  dslope = detrend(slope);
  warning('off','signal:findpeaks:largeMinPeakHeight')
  [valleys,valley_ix] = findpeaks(-dslope,'MinPeakHeight',minpeak,...
    'MinPeakDistance',10,'SortStr','descend');
  warning('on','signal:findpeaks:largeMinPeakHeight')  

  if isempty(valley_ix)
    k = -2; % No transition found
    return
  end 

  ixvalleys = valley_ix + start - 1;  % Index relative to start of trace
  validvalleys = valleys>valleys(1)*0.7;
  valleys = valleys(validvalleys);
  ixvalleys = ixvalleys(validvalleys);
  if sgn == -1
    too_soon = f(ixvalleys) <= -20;  % Do not accept refolding at high force
    ixvalleys(too_soon) = [];
    valleys(too_soon) = [];
  end

%% Examine up to ten of the deepest valleys

% The valleys are ordered from the deepest to more shallow.  We accept the deepest 
% that seems OK.
  px = polyfit(validrange,x(validrange),1);  % Fit straight line to x to combat noise
  dt = (t(end)-t(1))/(nf-1);  % Time step
  k = 0;
  ii = 1;

  %  *** Search among valleys to find the first valid(?) event ****
  while k<1 && ii <= min(numel(valleys),10)
    [k,force,dx,Fdot,ixfitb,pb,ixfita,pa,shift,noise] = ...
      check_event(ixvalleys(ii),f,start,stop,px,dt,sgn,algpar);
    ii = ii + 1;
  end
  if k<1
    return   % No valid event found
  end  
  dFdx = pb(1)/px(1);

%% Plot result
  if plotting
    figure;
    xfit = sgn*polyval(px,ixfitb(1):ixfita(end));
    xpos = sgn*polyval(px,k);        % x at transition
    subplot(411)
    plot(sgn*f);
    hold on;plot(k,sgn*f(k),'*r')
    ylabel('Force')
    title(s.file);
    subplot(412);
    plot(sgn*x)
    hold on;plot(k,sgn*x(k),'*r')
    xlabel('index') 
    ylabel xpos
    subplot(212)
    h = plot(sgn*x,sgn*f,xfit(1:numel(ixfitb)),sgn*polyval(pb,ixfitb),'k',...
      xfit(end-numel(ixfita)+1:end),sgn*polyval(pa,ixfita),'k');
    set(h(2:3),'linewidth',2)
    hold on;
    plot([xpos,xpos+sgn*dx],force*[1,1],'r','linewidth',2)
    hold off
    xlabel('xpos (nm)');
    ylabel('Force (pN)')
    title(sprintf('t = %.1fs, dx = %.1fnm',s.t(k),dx))
  end
end
%%
function [k,force,dx,Fdot,ixfitb,pb,ixfita,pa,shift,noise] = ...
    check_event(ixvalley,f,start,stop,px,dt,sgn,algpar)
  % Check if the event at ixvalley is a valid transition
  % Returns k = NaN if not valid
  force = NaN;
  dx = NaN;
  Fdot = NaN;
  pa = 0;
  pb = 0; 
  ixfita = [];
  ixfitb = [];
  shift = NaN;
  noise = NaN;

  nf = numel(f);
  % Find highest f value before ixvalley
  peakrange = max(ixvalley-10,1):ixvalley;
  [~,m] = max(f(peakrange));
  k = peakrange(1) - 1 + m;
  [~,m] = min(f( k+1:min(nf,k+12)));
  i2 = k+m;   % Lowest value after transition


  if k-start < algpar.minpointforfit || stop-i2 < algpar.minpointforfit
    % not enough points for fitting pa0 or pb0
    k = -3;
    return
  end      

  % Fit lines before and after (potential) transition
  fitL = round(nf/10);  % Fit length
  ixfitb = max(start,k-fitL):k;
  ixfita = k+3:min(k+fitL+3,stop);
  fb = f(ixfitb);
  fa = f(ixfita);
  pb = polyfit(ixfitb,fb,1);
  pa = polyfit(ixfita,fa,1);

  if algpar.multiple
    % With multiple unfoldings and refoldings on the same trace,
    % unfolded points before the unfolding and refolded points after
    % should be excluded from the fits
    % pa and pb are now considered tentative lines. Find the differences 
    % between measurement and tentative line
    dfb = fb-polyval(pb,ixfitb)'; 
    dfa = fa-polyval(pa,ixfita)';

    shift1 = polyval(pb,k)-polyval(pa,k);
    ok_b = dfb > -shift1*0.08;   % Eliminates points far from the line 
    ok_a = dfa < shift1*0.08;
    if sum(ok_a) < algpar.minpointforfit || sum(ok_b) <algpar.minpointforfit
      k = -4;  % To few relevant points for meaningful fit
      return
    end
    pb = polyfit(ixfitb(ok_b),fb(ok_b),1);
    pa = polyfit(ixfita(ok_a),fa(ok_a),1);
  else
    ok_b = logical(ones(size(ixfitb)));
  end

  if max(abs(pa(1)/pb(1)-1),abs(pb(1)/pa(1)-1)) > algpar.maxsloperatio
    k = -5;  % Slope before and after differ too much
    return
  end
  
  force = sgn*polyval(pb,k);   % Force at transition
  shift = polyval(pb,k)-polyval(pa,k);
  noise = 4*std(f(ixfitb(ok_b))-polyval(pb,ixfitb(ok_b))');
  if sgn > 0
    noiselimit = max(noise*0.7,0.8);
  else
    noiselimit = max(noise*0.4,0.5);
  end
  if shift < noiselimit
    k = -6;  % Too small shift
    return
  end

  % polyval(pa,recovery_ix) == force
  recovery_ix = (sgn*force-pa(2))/pa(1); 
  dx = px(1)*(recovery_ix-k);
  Fdot = pb(1)/dt;
  if abs(dx) < 5 | abs(dx) > 30
    k = -7;   % unrealistic dx
  end
end