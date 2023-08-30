  close all
if ~exist('lidar_sm')
  load ./data/summary_lidar_science_small.mat
end
if ~exist('allresults')
  load allresults
end
for j  = 1:length(allresults)
  for i = 1:length(allresults(j).results.hydro)
    [junk ind ]=min(abs([lidar_sm.date]-allresults(j).g.date_bc(i)));
    allresults(j).lid_r2p(i) = lidar_sm(ind).r2p;
    allresults(j).lid_Hs(i,:) = lidar_sm(ind).wg_Hs;
  end
end
x  = lidar_sm(ind).wg_x;
xcshore = allresults(1).results(1).hydro(1).x;
x2 = allresults(1).g.x_offset-allresults(1).results(1).hydro(1).x;
for i  = 1:length(lidar_sm(ind).wg_Hs);
  [ junk wginds(i)] = min(abs(x(i)-x2));
end

fs = 14;

for j  = 1:length(allresults)
  figure(j);clf
  hh(1)=plot(allresults(j).g.date_bc,allresults(j).lid_r2p,'rs','markerfacecolor','k');hold all
  hh(2)=plot(allresults(j).g.date_bc,[allresults(j).results.hydro.runup_2_percent],'b-','linewidth',2);hold all
  hh(3)=plot(allresults(j).g.date_bc,[allresults(j).g.wl],'k-','linewidth',2);hold all
  
  title('$R_{2\%}$ Runup','fontsize',fs,'interpreter','latex')
  ylabel('$R_{2p}[m]$','fontname','times','fontsize',fs,'interpreter','latex')
  xlabel('date','fontname','times','fontsize',fs,'interpreter','latex')
  datetick
  legend(hh,'Data','CSHORE','Water level at SB')
    set(gcf,'position',[600 500 900 400])
  print('-dpng', ['runup_',num2str(j),'.png']);
end 

for j  = 1:length(allresults)
  Hs = [allresults(j).results.hydro.Hrms]*sqrt(2);
  figure(j+length(allresults));clf;orient tall;
  for k = length(x):-1:1
    subplot(length(x),1,length(x)-k+1);
    plot(allresults(j).g.date_bc,Hs(wginds(k),:),'b-','linewidth',2,'markerfacecolor','k');hold all
    plot(allresults(j).g.date_bc',allresults(j).lid_Hs(:,k),'rs','markerfacecolor','k');hold all
    a = axis;
    text(a(1)+.01*(a(2)-a(1)),a(4)-.2*(a(4)-a(3)),{['$x_{frf}$ = ',num2str(x(k))]; ...
                  ['$x_{cshore}$ = ',num2str(xcshore(wginds(k)))]; },'interpreter','latex') 
    datetick
    ylabel('$H_{s}[m]$','fontname','times','fontsize',fs,'interpreter','latex')
    if k==length(x);title('$H_{s}$ Wave Height','fontsize',fs,'interpreter','latex');end
    if k~=1;set(gca,'xticklabels',[]);end
    if k==1;xlabel('date','fontname','times','fontsize',fs,'interpreter','latex');end

  end
  
    print('-dpng', ['Hs_',num2str(j),'.png']);
  %  hh(2)=plot(allresults(j).g.date_bc,[allresults(j).results.hydro.runup_2_percent],'b-','linewidth',2);hold all
  %hh(3)=plot(allresults(j).g.date_bc,[allresults(j).g.wl],'k-','linewidth',2);hold all
  

  %legend(hh,'Data','CSHORE','Water level at SB')
end 


for j  = 1:length(allresults)
  figure
  x = allresults(j).g.x_offset-allresults(j).results.morpho.x;
  x2  = lidar_sm(ind).wg_x;
  zb = allresults(j).results.morpho.zb;
  plot(x,zeros(size(x)),'linewidth',3);hold all
  fill([x;x(end);x(1)],[zb; -20; -20],[.9 .9 .7])
  plot([x2;x2],[ones(size(x2));-ones(size(x2))],'k')
  axis([70 600 -6 2.9])
  set(gcf,'position',[600 500 900 400])
  title(['FRF Domain, ',datestr(allresults(j).g.begin_date)],'fontsize',fs,'interpreter','latex')
  ylabel('$z [m]$','fontname','times','fontsize',fs,'interpreter','latex')
  xlabel('x [m]','fontname','times','fontsize',fs,'interpreter','latex')
  text(100,2,{'Inner Surf';'Wave Gauges'},'fontname','times','fontsize',fs,'interpreter','latex')
    print('-dpng', ['domain_',num2str(j),'.png']);
end

  
  





