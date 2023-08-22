function [k,force,dx,Fdot,shift,pullingspeed] = analyse_stretch(s,plotting)
% Find transition point and distance in a unfolding or refolding stretch
% Inputs;:
%  s:  input struct with fields:
%    x: extent (trap position)
%    f: force 
%    t: time
%  plotting: plot f(x) with linear fits and dx
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
%  -movingslope by John D'Errico
%    https://www.mathworks.com/matlabcentral/fileexchange/16997-movingslope
%    MATLAB Central File Exchange. Retrieved October 2, 2022.

% Version 2.12
% Author: Are Mjaavatten (mjaavatt@gmail.com)
% Date:   2023-08-19

% Change history
% v 2.8:  Extra output: shift
%         No rounding of recovery_ix.  Avoids spurious discretization of dx
% v 2.9   Restructured to reflect analyse_stretch_demo.mlx
% v 2.10  Added ouput: pullingspeed
% v 2.11  Corrected typos.
% v 2.12  algo_par may override standard algorithm paramters

% First step in plan to collect algorithm tuning parameters in one place
  legpar.minforcerange = 10;
  legpar.minpointforfit = 3;
  legpar.maxsloperatio = 0.5;

  if nargin < 2
    plotting = false;
  end
%% Tuning parameters
  % Transitions are unlikely near the start of a stretch, so we exclude
  % potential transitions if abs(f-fstart) < legpar.dfstart
  legpar.dfstart = 4;  
  % At the end of a REFOLD stretch noise makes detections difficult:
  legpar.dfend = 1;  % Make sure that abs(f-fend) > legpar.dfend at transition.
  % No end restrictions for unfold stretches

  % Low slope parts at start or end of a stretch may confuse the detection 
  % algorithm, so eliminate start and end parts where the slope
  % is less that legpar.slopefrac*median(slope):
  legpar.slopefrac = 0.7;

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

% Make sure x>0 and f both increase or decrease together
  s.x = s.x - min(s.x);
  
% Change sign of x and f for refolding stretches
% This lets us use the same algorithm for bot unfolding and refolding
  sgn = sign(s.f(end) - s.f(1));
  f = sgn*s.f;
  x = sgn*s.x;
  t = s.t;

% Do not accept very short stretches
  df = f(end)-f(1);  % range of forces
  if(df<legpar.minforcerange)  % Skip if force amplitude is small (probably noise only)
    k = -1;
    return
  end

% Eliminate bad parts at start and end
  slope = movingslope(f,slopelength);
  medslope = median(slope);
  startslope = find(slope>legpar.slopefrac*medslope,1);  % End of inital flat part 
  startmin = find(abs(f(1)-f)>legpar.dfstart,1);
  start = max([1,startslope,startmin]);
  stopslope = find(slope>legpar.slopefrac*median(slope),1,'last'); 
  stopmin = find(f(nf) - f > legpar.dfend,1,'last');
  stop = min([stopmin,stopslope,nf]);
  validrange = start:stop;
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

  ixvalleys = valley_ix + start - 1;  % Index relative to start of stretch
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
    [k,force,dx,Fdot,ixfitb,pb,ixfita,pa,shift] = ...
      check_event(ixvalleys(ii),f,start,stop,px,dt,sgn,legpar);
    ii = ii + 1;
  end
  if k<1
    return   % No valid event found
  end  

%% Plot result
  if plotting
    xfit = sgn*polyval(px,ixfitb(1):ixfita(end));
    xpos = sgn*polyval(px,k);        % x at transition
    h = plot(sgn*x,sgn*f,xfit(1:numel(ixfitb)),sgn*polyval(pb,ixfitb),'k',...
      xfit(end-numel(ixfita)+1:end),sgn*polyval(pa,ixfita),'k');
    set(h(2:3),'linewidth',2)
    hold on;
    plot([xpos,xpos+sgn*dx],force*[1,1],'r','linewidth',2)
    hold off
    xlabel('xpos (nm)');
    ylabel('Force (pN)')
    title(s.file); 
  end
end

function [k,force,dx,Fdot,ixfitb,pb,ixfita,pa,shift] = ...
    check_event(ixvalley,f,start,stop,px,dt,sgn,legpar)
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

  nf = numel(f);
  % Find highest f value before ixvalley
  peakrange = max(ixvalley-10,1):ixvalley;
  [~,m] = max(f(peakrange));
  k = peakrange(1) - 1 + m;
  [~,m] = min(f( k+1:min(nf,k+12)));
  i2 = k+m;   % Lowest value after transition


  if k-start < legpar.minpointforfit || stop-i2 < legpar.minpointforfit
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

  % ******* Trick to handle multiple unfoldins and refoldings  *****
  pb1 = polyfit(ixfitb,fb,1);     % Tentative line before
  pa1 = polyfit(ixfita,fa,1);
  dfb = fb-polyval(pb1,ixfitb)';  % Difference between measurement and tentative line
  dfa = fa-polyval(pa1,ixfita)';

  %E.g.: For unfolding, fit to only unfolded points before the potential
  % transition and only folded points after
  shift1 = polyval(pb1,k)-polyval(pa1,k);
  ok_b = dfb > -shift1*0.08;   % Eliminates points far from the line 
  ok_a = dfa < shift1*0.08;

  % ******  Undo the trick ******
  % This should be handled more elegantly in a later version !!!!
  ok_a = ok_a | 1;
  ok_b = ok_b | 1;


  if sum(ok_a) < legpar.minpointforfit || sum(ok_b) <legpar.minpointforfit
    k = -4;  % To few relevant points for meaningful fit
    return
  end

  pb = polyfit(ixfitb(ok_b),fb(ok_b),1);
  pa = polyfit(ixfita(ok_a),fa(ok_a),1);
  if max(abs(pa(1)/pb(1)-1),abs(pb(1)/pa(1)-1)) > legpar.maxsloperatio
    k = -5;  % Slope before and after differ too much
    return
  end
  
  force = sgn*polyval(pb,k);   % Force at transition
  shift = polyval(pb,k)-polyval(pa,k);
  noise = 4*std(f(ixfitb(ok_b))-polyval(pb,ixfitb(ok_b))');
  if sgn == 1
    % noiselimit = max(noise,0.6+0.4*noise);  % higher limit if noise is low
    noiselimit = max(noise,0.35+0.65*noise);  % Ad-hoc override
  else
    % noiselimit = max(noise,0.35+0.65*noise);  % shift is lower for refodling
    noiselimit = max(noise,0.25+0.65*noise);
  end
 
  if shift < noiselimit
    k = -6;  % Too small shift
    return
  end

  % polyval(pa,recovery_ix) == force
  recovery_ix = (sgn*force-pa(2))/pa(1); 
  dx = px(1)*(recovery_ix-k);
  Fdot = pb(1)/dt;
  if abs(dx) < 10
    k = -7;   % Too big dx
  end
end