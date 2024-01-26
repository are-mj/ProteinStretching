function [DG,DGstd] = fit_Crooks(Tun,Tre,ucases,rcases,Clusters,plotting)
% Calculate DG from unfolding and refolding parameters dx and log10k0 
% Unfolding for high and low force clusters, Refolding for all data.
% DG equals the work at the point where unfolding and refolding work are
% equal, according to the Crooks fluctuation theorem
% Output:
%   DG:    Crooks DeltaG
%   DGstd: Standard deviation of DG calculated by MonteCarlo simulation
%          Time-consuming. Calculation skipped if nargout < 3

  if nargin < 5
    plotting = 0;
  end
  P = 0.65;L0 = 29.28;  % WLC parameters

  Tmean = mean(Tun.Temperature(ucases.selected));
  T = Tmean + 273.15;

  % Refold:
  f_refold = Tre.Force(rcases.selected);
  deltax_refold = Tre.Deltax(rcases.selected);
  Ws = stretchwork(f_refold,deltax_refold,P,T,L0); % Should we subtract Ws here?
  Wr = f_refold.*deltax_refold - Ws;
  Wr_kcal = Wr*0.1439;  % Convert from nm*pN per molecule to kcal/mol
  pdr = fitdist(Wr_kcal,'normal');   % Normal distribution

  % Unfold:
  f_unfold = Tun.Force(ucases.selected);
  deltax_unfold = Tun.Deltax(ucases.selected);

  tbl = Tun;
  selected = ucases.selected;
  n = sum(Clusters);
  nclusters = find(n>8);    % Clusters with more that 8 elements  

  % Tmean = mean(tbl.Temperature(selected),'omitnan');
  % Fdot = mean(tbl.Fdot(selected));

  ncl = numel(nclusters);   % Number of clusters
  fu = cell(ncl,1);
  dxu = cell(ncl,1);
  Wu = cell(ncl,1);
  pdu = cell(ncl,1);
  DG = zeros(2,1);
  DGstd = zeros(2,1);
  if plotting
    hh = [];
    nexttile;
    h = plot(pdr);hold on
    hh = [hh,h];
    colors = [0 0 1;0.5 0.2 0.5;0.3 0.6 0.6];
    h(1).Color = colors(1,:);
    % legtext = {'Refolding data','Normal distribution'};
    legtext = "Refolding data";
  end
  for cl = find(n>8)  % Nonzero parameter columns
    fu{cl} = f_unfold(Clusters(ucases.selected,cl));
    dxu{cl} = deltax_unfold(Clusters(ucases.selected,cl));
    Ws = stretchwork(fu{cl},dxu{cl},P,T,L0);
    Wu{cl} = (fu{cl}.*dxu{cl}-Ws)*0.1439;        % Net work
    pdu{cl} = fitdist(Wu{cl},'normal');
    DG(cl) = match(pdr,pdu{cl});
    if nargout > 2
      % Perform this time-consuming calculation only if needed:
      DGstd(cl) = Crooks_std(pdu{cl},pdr);
    end
    if plotting
      h = plot(pdu{cl});
      hh = [hh,h];
      h(1).Color = colors(cl+1,:);
      if cl == 1
        s = 'Low force unfolding';
      else
        s = 'High force unfolding';
      end
      % legtext = [legtext,{s,'Normal distribution'}];
      legtext = [legtext;s];
      title(ucases.text)
      xlabel('');ylabel(''); % remove fitdist's standard labels
    end
  end
  if plotting
    set(gca,'XScale','log')
    xlim([2,200]);
    ylim([0,0.3]);
    set(gca,'xtick',[3,10,30,100])
    legend(hh(2:2:6),legtext);
  end
end

function pd = pdfun(pdobj,x)
  mu = pdobj.mu;
  sigma = pdobj.sigma;
  pd = exp(-0.5*((x-mu)/sigma).^2)/sigma/sqrt(2*pi);
end

function DG = match(pd1,pd2)
% Find the point where two normal distributions are equal
  fun = @(x) pdfun(pd1,x)-pdfun(pd2,x);
  % [DG,~,exitflag] = fzero(fun,(pd1.mu+pd2.mu)/2);
  x0 = 0;
  while fun(x0)<0  &&  x0 < pd1.mu+pd2.mu
    x0 = x0 + 5;
  end
  [DG,~,exitflag] = fzero(fun,[x0,pd1.mu+pd2.mu]);
  if exitflag < 0
    warning('fzero problem')
  end
end

function DGstd = Crooks_std(pdu,pdr)
% Monte Carlo calculation of the std of the intersection of pdr and pdu
  
  pdu_ci = pdu.paramci;
  pdu_mu_std = pdu.mu - pdu_ci(1,1);
  pdu_sigma_std = pdu.sigma - pdu_ci(1,2);

  pdr_ci = pdr.paramci;
  pdr_mu_std = pdr.mu - pdr_ci(1,1);
  pdr_sigma_std = pdr.sigma - pdr_ci(1,2);

  n = 10000;
  % n = 100;  % Faster, but less accurate

  DG = zeros(n,1);
  for i = 1:n
    try  % normrnd may occasionally return negative values
         % This will crash makedist
      umu = normrnd(pdu.mu,pdu_mu_std);
      usigma = normrnd(pdu.sigma,pdu_sigma_std);
      rand_pdu = makedist('normal','mu',umu,'sigma',usigma);
      % plot(rand_pdu);
  
      rmu = normrnd(pdr.mu,pdr_mu_std);
      rsigma = normrnd(pdr.sigma,pdr_sigma_std);
      rand_pdr = makedist('normal','mu',rmu,'sigma',rsigma); 
      % plot(rand_pdr);
      DG(i) = match(rand_pdr,rand_pdu);
    catch
      DG(i) = NaN;
    end
  end
  DG(isnan(DG)) = [];
  DGstd = std(DG);
end
