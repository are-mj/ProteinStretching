function [dx1,dx2,t] = unfold_distances(file,plotting)
  % dx1: unfold distance, dx2 = distance to refold at unfold pressure

  if nargin < 2
    plotting = 0;
  end
  [Tu,Tr,~,f0,x0,~,peaks] = identify_events(file);
  
  x0 = -x0 + max(x0);  % Make sure x increases during stretching 
  r = valid_ranges(file);
  if isempty(r)
    r = [1,numel(f0)];
  end
  dx1 = [];
  dx2 = [];
  t =[];
  for part = 1:size(r,1) 
    f = f0(r(part,1):r(part,2));
    x = x0(r(part,1):r(part,2));
    [lopos,hipos,Flo,Fhi] = necklace(peaks(part).lopos,peaks(part).hipos);
    peaks(part).lopos = lopos(:);
    peaks(part).hipos = hipos(:);
    peaks(part).hival(Fhi) = [];
    peaks(part).loval(Flo) = [];
    bad = peaks(part).hival - peaks(part).loval(1:end-1) < 10;  % Eliminate peaks with small amplitude
    peaks(part).hival(bad) = [];
    peaks(part).hipos(bad) = [];
    peaks(part).loval(bad) = [];
    peaks(part).lopos(bad) = [];
  
    for ii = 1:numel(peaks(part).hipos)
      % Indices for valley-peak-valley sequence 
      i0 = peaks(part).lopos(ii);    % Valley 1
      i1 = peaks(part).hipos(ii);    % Peak
      i2 = peaks(part).lopos(ii+1);  % Valley2
   
      Tu_row = find(Tu.Lineno>i0+r(part,1)-1 & Tu.Lineno<i1+r(part,1)-1);
      Tr_row = find(Tr.Lineno>i1+r(part,1)-1 & Tr.Lineno<i2+r(part,1)-1);   
      if isempty(Tu_row) || isempty(Tr_row)
        continue
      end
      k1 = Tu.Lineno(Tu_row) - r(part,1); 
      x1 = x(k1);
      k2 = Tr.Lineno(Tr_row) - r(part,1);     
      f2 = f(i1:k2);
      x2 = x(i1:k2);
      
      % Moving mean filter
      nfilter = round((i2-i1)/25); 
      b = 1/nfilter*ones(nfilter,1);  % Numerator filter coefficients
      a = 1;  % Denominator
      x2filt = filtfilt(b,a,x2);
      f2filt = filtfilt(b,a,f2);
      stop = numel(f2filt);
      start = find(diff(f2filt)>=0,1,'last')+2;
      if isempty(start) | start/numel(f2filt) > 0.7
        start = 1;
        stp = find(diff(f2filt)>=0,1)-2;
        if ~isempty(stp)
          stop = stp;
        end
      end
      x2 = interp1(f2filt(start:stop),x2filt(start:stop),Tu.Force(Tu_row));
      dx1 = [dx1,Tu.Deltax(Tu_row)];
      dx2 = [dx2,x2-x1];
      t = [t,Tu.Time(Tu_row)];
    end
  end
  if plotting
    figure;
    set(gcf,'Position',[2000,500,400,280])
    semilogy(t,dx1,t,dx2,t);
    ylim([0.5,100])
    title(file);
    xlabel Time(s)
    legend('\Deltax1','\Deltax2','location','best')  
    title(file)
  end
end
