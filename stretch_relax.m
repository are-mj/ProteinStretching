% stretch_relax.m
% Utility for extracting a stretch-relax curve from a plot with t on x axis
% and title equal to filename (e.g.'07182022/cAB.txt')
  

  show_mean_xx0 = 0;
  try
    file = get(gca).Title.String;
    [shortname,bead,date,t0,f0,xx0,T] = read_experiment_file(file);
    x0 = mean(xx0,2);
    x0 = -x0 + max(x0);  % Make sure x increases during stretching 
    xx0 = -xx0+ max(xx0);
  catch
    error('The title of the current figure must be a valid file name')
  end
  uiwait(msgbox('Select a time'));
  [tt,~] = ginput(1);  
  ii = round(interp1(t0,1:numel(t0),tt));
    % Estimate length of cycle:
  fm = median(f0);
  fmhi = median(f0(f0>fm));
  fmlo = max(10,median(f0(f0<fm)));
  next_hi = find(f0(ii:end)>fmhi,1)+ii+1;
  next_lo = find(f0(ii:end)<fmlo,1)+ii+1;
  last_hi = find(f0(1:ii)>fmhi,1,'last');
  last_lo = find(f0(1:ii)<fmlo,1,'last');

  deltai = (max(next_hi,next_lo)-min(last_hi,last_lo))*3;
  testrange = max(1,ii-deltai):min(ii+deltai,numel(t0));
  figure; plot(testrange,f0(testrange),ii,f0(ii),'*r')
  uiwait(msgbox('Select start of stretch_relax pair'));
  [i0,~] = ginput(1);
  uiwait(msgbox('Select end of stretch_relax pair'));
  [i2,~] = ginput(1); 
  i0 = round(i0);
  i2 = round(i2);
  [~,j1] = max(f0(i0:i2));
  i1 = i0+j1-1;
  delta_ix = 0;  % Delay in recording x (detat_t = delta_ix/200)
  figure;
  if show_mean_xx0
    plot(x0((i0:i1)+delta_ix),f0(i0:i1),x0((i1:i2)+delta_ix),f0(i1:i2))
  else
    plot(xx0((i0:i1)+delta_ix,1),f0(i0:i1),xx0((i1:i2)+delta_ix,1),f0(i1:i2))
  end
  xlabel x;ylabel f
  title(sprintf('%s.   t = %.1f - %.1fs',file,t0(i0),t0(i2)))
  xlabel x;ylabel f
  legend('stretching','relaxing')

