function [ok,Cluster1,Cluster2,par] = no_outliers(Tun,par,plotting)
% Remove outliers as defined by dbscan, keeping only Cluster1 and Cluster2
% ok:       logical array, true for non-outlier rows in Tun
% Cluster1: logical array, true for Cluster1 rows in Tun
% Cluster2: logical array, true for Cluster2 rows in Tun 
%
%  (ok = Cluster | Cluster2,  outliers = ~ok)

  if nargin < 2 || isempty(par)
    par.epsilon = 5.5;
    par.minpts = 200;
    par.scaling = 3;
  end


  X = [Tun.Deltax*par.scaling,Tun.Force];  % Scaled
  labels = dbscan(X,par.epsilon,par.minpts);
  Cluster1 = labels == 1;
  Cluster2 = labels == 2;
  ok = labels>0;
  if nargin>2
    if plotting
      figure;
      plot(Tun.Deltax(Cluster1),Tun.Force(Cluster1),'.b', ...
        Tun.Deltax(Cluster2),Tun.Force(Cluster2),'.r');
      hold on;
      plot(Tun.Deltax(~ok),Tun.Force(~ok),'.','color',0.65*[1 1 1])
      hold off
    end
  end
end