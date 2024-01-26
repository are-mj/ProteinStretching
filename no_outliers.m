function [ok,Cluster1,Cluster2,par] = no_outliers(Tun)
% Remove outliers as defined by dbscan, keeping only Cluster1 and Cluster2
% ok:       logical array, true for non-outlier rows in Tun
% Cluster1: logical array, true for Cluster1 rows in Tun
% Cluster2: logical array, true for Cluster2 rows in Tun 
%
%  (ok = Cluster | Cluster2,  outliers = ~ok)

  par.epsilon = 5.5;
  par.minpts = 200;
  par.scaling = 3;
  X = [Tun.Deltax*par.scaling,Tun.Force];  % Scaled
  labels = dbscan(X,par.epsilon,par.minpts);
  Cluster1 = labels == 1;
  Cluster2 = labels == 2;
  ok = labels>0;
end