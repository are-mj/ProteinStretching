function s = display_trace(tbl,xcolumn,ycolumn,fig)
% Display the traces defined by the data tips in the current figure
%  tbl: table from analyse_many
%  xcolumn: Column in tbl on x axis. Default: Deltax
%  ycolumn: Column in tbl on y axis. Default: Force
%  Example: Create figure by
%   Tfast = select(Tun,'Pullingspeed>200');
%   xx = linspace(10,25); figure; 
%   plot(Tfast.Deltax,Tfast.Force,'.',xx,wlc(xx,.65,290,29.28),'k')
%  Create data tips
%  display_trace(Tfast)
  
  if nargin < 4
    fig = gcf;  % Use current figure if 
  end
  if nargin < 3
    ycolumn = 'Force';
  end
  if nargin < 2
    xcolumn = 'Deltax';
  end
  [x,y] = read_data_tip(fig);
  s = [];
  for i = 1:length(x)
    T = tbl(tbl.(xcolumn) == x(i) & tbl.(ycolumn) == y(i),:);
    s = [s;extract_trace(T.Filename,T.Time)];
    figure;
    subplot(211);
    plot(s(i).x,s(i).f);
    subplot(212);
    analyse_trace(s(i),1);
  end
end

function [x,y] = read_data_tip(fig)
% Return x and y values for all data tips in current figure
  if nargin < 1
    fig = gcf;
  end
  dcm = datacursormode(fig);
  cursorstructs = getCursorInfo(dcm);
  data = [cursorstructs.Position]';
  x = data(1:2:end);
  y = data(2:2:end);
end