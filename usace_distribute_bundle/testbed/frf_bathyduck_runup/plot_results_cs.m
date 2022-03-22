close all
if ~exist('lidar')
  load ./data/summary_lidar_science.mat
end
if ~exist('allresults')
  load allresults
end
for j  = 1:length(allresults)
  for i = 1:length(allresults(j).results.hydro)
    [junk ind ]=min(abs([lidar.date]-allresults(j).g.date_bc(i)));
    allresults(j).lid_r2p(i) = lidar(ind).r2p;
    allresults(j).lid_Hs(i,:) = lidar(ind).wg_Hs;
  end
end
x  = lidar(ind).wg_x;
xcshore = allresults(1).results(1).hydro(1).x;
x2 = allresults(1).g.x_offset-allresults(1).results(1).hydro(1).x;
for i  = 1:length(lidar(ind).wg_Hs);
  [ junk wginds(i)] = min(abs(x(i)-x2));
end

fs = 24;

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
end 

for j  = 1:length(allresults)
  Hs = [allresults(j).results.hydro.Hrms]*sqrt(2);
  figure(j+length(allresults));clf;orient tall;
  for k = length(x):-1:1
    subplot(length(x),1,length(x)-k+1);
    plot(allresults(j).g.date_bc,Hs(wginds(k),:),'b-','markerfacecolor','k');hold all
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
  
  
  %  hh(2)=plot(allresults(j).g.date_bc,[allresults(j).results.hydro.runup_2_percent],'b-','linewidth',2);hold all
  %hh(3)=plot(allresults(j).g.date_bc,[allresults(j).g.wl],'k-','linewidth',2);hold all
  

  %legend(hh,'Data','CSHORE','Water level at SB')
end 




