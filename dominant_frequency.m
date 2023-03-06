function f = dominant_frequency(t,y)
% Find the dominant frequancy in a time series
%  t: time vector
%  y: value vector
%  f: Dominant frequency
% Based on Star Strider's accepted answer in:
%   https://se.mathworks.com/matlabcentral/answers/460152

  offset = 10;  % Avoid spurious low frequaencies
  Ts = mean(diff(t));
  Fs = 1/Ts;
  Fn = Fs/2;
  L = numel(t);
  ym = y-mean(y);                    % Subtract Mean (Mean = 0 Hz)
  FCl = fft(ym)/L;
  Fv = linspace(0, 1, fix(L/2)+1)*Fn;
  Iv = offset:numel(Fv);
  [v,locs] = findpeaks(abs(FCl(Iv))*2, 'MinPeakHeight',0.01);
  [~,m] = max(v);
%   figure
%   plot(Fv, abs(FCl(Iv))*2)
%   grid
%   xlim([0 0.1])
  f = Fv(locs(m)+offset-1);
end
