function [split,outbell,outdudko] = fit_dual(tbl,casestruct)
% Fit parameters to data, assuming two force clusters
% Input:
%   tbl:  Matlab table Tun or Tre from analyse_many (or analys_table)
%   casestruct: struct specifying the the case to study. Specified in
%     run_dual_model.  casestruct has the following fields: 
%        selected:  logical column array with as many rows as tbl
%                     true for selected rows
%        speedtest: string describing the pulling speed for the selection
%        temptext:  string describing the temperatures selected
%

  % Allocate variables
  unfold = tbl.Fdot(1) > 0;
  if unfold
    Fplot = linspace(5,55)';
    dF = 1;   % Force histogram spacing
    split0 = 22.5;
    thetabell0 = [1;-2];
  else
    Fplot = linspace(3,18)';
    dF = 0.5;   % Force histogram spacing
    split0 = 8;    
    thetabell0 = [4;4];
  end
  nplot = numel(Fplot);
  pdbell = NaN(nplot,2);
  pddudko = NaN(nplot,2);
  pdplot = NaN(nplot,2);
  mark = zeros(2);
  outbell = zeros(2,6);
  outdudko = zeros(2,8);

  speedtext = casestruct.speedtext;
  temptext = casestruct.temptext;
  selected = casestruct.selected;
  Tmean = mean(tbl.Temperature(selected));
  Fdot = mean(tbl.Fdot(selected));

  force = tbl.Force(selected);

  fun = @(split) residual_bell(split,force,Tmean,Fdot,thetabell0);
  split  = fminsearch(fun,split0);

  cluster1 = force<=split;
  cluster2 = force>split;
 
  % Dudko model inputs:
  par.nu = 0.5;
  par.model = 'DHS';
  thetaD0 = [50;5;-4];
  
  [pd{1},F{1},n(1)] = probdens(force(cluster1),dF);
  [pd{2},F{2},n(2)] = probdens(force(cluster2),dF);

  if all(n>4)  % No empty cluster
    for cl = 1:2     
      if unfold
      [bellth(:,cl),bellstd(:,cl),bellnorm(1,cl)] = fit_Bell_unfold(F{cl},pd{cl},Tmean,Fdot,thetabell0);
      pdbell(:,cl) = Bell_unfold_probability(bellth(:,cl),Fplot,Tmean,Fdot);          
        [dudkoth(:,cl),dudkostd(:,cl),dudkonorm(1,cl)] ...
          = fit_Dudko_unfold(F{cl},pd{cl},Tmean,Fdot,thetaD0,par);
        thetacalc = [dudkoth(1,cl);par.nu*dudkoth(2,cl)/dudkoth(1,cl);dudkoth(3,cl)];
        pddudko(:,cl) = Dudko_unfold_probability(thetacalc,Fplot,Tmean,Fdot,par);
        bellbest(cl) = bellnorm(cl)./dudkonorm(cl) < 0.96 | any(dudkostd(:,cl)>2*max(abs(dudkoth(:,cl)),1));
        pdplot(:,cl) = bellbest(cl)*pdbell(:,cl) + ~bellbest(cl)*pddudko(:,cl);  
        mark(cl,:) = bellbest(cl).*double('* ') + ~bellbest(cl).*double(' *'); 
      else
        [bellth(:,cl),bellstd(:,cl),bellnorm(1,cl)] = fit_Bell_refold(F{cl},pd{cl},Tmean,Fdot,thetabell0);
        pdbell(:,cl) = Bell_refold_probability(bellth(:,cl),Fplot,Tmean,Fdot);          
        pdplot = pdbell;
      end
      print_results;
    end
    plot_results;
  else  % Only one cluster
    cl = find(n>4);
    if unfold
      [bellth,bellstd,bellnorm] = fit_Bell_unfold(F{1},pd{1},Tmean,Fdot,thetabell0);
      pdbell(:,cl) = Bell_unfold_probability(bellth,Fplot,Tmean,Fdot);
      [dudkoth,dudkostd,dudkonorm] = fit_Dudko_unfold(F{1},pd{1},Tmean,Fdot,thetaD0,par);
      thetacalc = [dudkoth(1);par.nu*dudkoth(2)/dudkoth(1);dudkoth(3)];
      pddudko(:,cl) = Dudko_unfold_probability(thetacalc,Fplot,Tmean,Fdot,par);
     
      bellbest = bellnorm./dudkonorm < 0.96 | any(dudkostd>2*max(abs(dudkoth),1))';
      pdplot = bellbest*pdbell + ~bellbest*pddudko;  
      mark = bellbest*double('* ') + ~bellbest*double(' *'); 
    else
      [bellth(:,cl),bellstd(:,cl),bellnorm(1,cl)] = fit_Bell_refold(F{cl},pd{cl},Tmean,Fdot,thetabell0);
      pdbell(:,cl) = Bell_refold_probability(bellth(:,cl),Fplot,Tmean,Fdot);         
      pdplot = pdbell;
    end
    print_results;
    plot_results
  end

  function print_results
    out = [bellth(:,cl)';bellstd(:,cl)'];
    outbell(cl,:) = [out(:)',n(cl),bellnorm(cl)];
    fprintf('%7s %7s %7d  Bell',speedtext,temptext,cl)
    fprintf('%1s                 %5.2f ±%4.1f %6.2f ±%4.1f    %4d %10.4f\n',char(mark(cl,1)),outbell(cl,:))
    if unfold
      out = [dudkoth(:,cl)';dudkostd(:,cl)'];
      outdudko(cl,:) = [out(:)',n(cl),dudkonorm(cl)];
      fprintf('%7s %7s %7d Dudko',speedtext,temptext,cl)
      fprintf('%1s   %6.2f ±%5.1f %5.2f ±%4.1f %6.2f ±%4.1f    %4d %10.4f\n',char(mark(cl,2)),outdudko(cl,:))    
    end
  end

  function plot_results
    figure;
    [pd,F] = probdens(force,dF);
    bar(F,pd,'FaceColor',[0.2,.7,.9],'BarWidth',1)
    hold on
    pdplot(isnan(pdplot)) = 0;
    w = n/sum(n);
    plot(Fplot,pdplot(:,1)*w(1)+pdplot(:,2)*w(2),'r','linewidth',2);
    plot(Fplot,pdbell(:,1)*w(1),'b',Fplot,pdbell(:,2)*w(2),'b');
    plot(Fplot,pddudko(:,1)*w(1),'.k',Fplot,pddudko(:,2)*w(2),'.k');
    title(strcat(speedtext," pulling. ",temptext," temperature","."))
    xlabel('Force (pN)')
    ylims= ylim;

    plot(split*[1,1],ylims,'--k')
    h = flip(get(gca,'children'));
    if unfold
      legend(h([1,2,3,5,7]),'Experiment','Dual model','Bell','Dudko','Force split')
    else
      legend(h([1,2,3,7]),'Experiment','Dual model','Clusters','Force split')
    end
    set(gca,"XTick",0:5:55);
    drawnow;
    hold off      
  end
end

function res = residual_bell(split,force,Tmean,Fdot,theta0)
  % theta0 = [1;-2];
  dF = 1;
  if Fdot < 0  % Refolding
    % theta0 = [5;5];
    dF = 0.5;
  end
  cluster{1} = force <= split;
  cluster{2} = force > split;
  resnorm = zeros(2,1);
  for cl = 1:2
    [pd,F,n(cl)] = probdens(force(cluster{cl}),dF);
    if n(cl) < 3
      resnorm(cl) = 10*n(cl);  % Trick to empty cluster with few points
    else
      if Fdot > 0
        [~,~,resnorm(cl)] = fit_Bell_unfold(F,pd,Tmean,Fdot,theta0);
      else
        [~,~,resnorm(cl)] = fit_Bell_refold(F,pd,Tmean,Fdot,theta0);
      end
    end
    
  end
  res = sum(resnorm);
end
