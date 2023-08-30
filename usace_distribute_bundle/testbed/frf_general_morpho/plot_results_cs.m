iprint = 1;
load allresults
fs = 18;
for j = 1:length(allresults)
  g=allresults(j).g;
  results=allresults(j).results;

  figure(j);clf
  x = [g.x_offset-g.x_frf]';
  hh(1) = plot(x,g.zb_init,'k','linewidth',4);hold all
  fill([x x(end) x(1)], [g.zb_init g.zb_init(end)-10 g.zb_init(1)-10],[.8 .8 .5])
  hh(2) = plot(x,g.zb_fin,'r','linewidth',4);
  hh(3) = plot(results.morpho(end).x,results.morpho(end).zb,'b--','linewidth',4);
  plot(results.hydro(1).x,0*results.hydro(1).x,'k')
  % plot(results_new.hydro(1).x,500*mean([results_new.hydro.dbsta],2))
  % plot(results_new.hydro(1).x,500*results_new.hydro(325).dbsta)
  % plot(results_new.sed(1).x_cross,5000*mean([results_new.sed.vs],2),'r--')
  axis([50 650 -6.5 6])
  set(gcf,'Position',[600 501 1000 400])
  hlabs{1} = datestr(g.begin_date);
  hlabs{2} = datestr(g.end_date);
  hlabs{3} = ['CSHORE'];
  legend(hh,hlabs,'location','northwest')
  title(['FRF Transect, y = ',num2str(g.y)],'interpreter','latex','fontsize',fs)
  ylabel('$z [m]$','interpreter','latex','fontsize',fs)
  xlabel('$x [m]$','interpreter','latex','fontsize',fs)

  if iprint;print('-dpng','-r300',['morpho',datestr(allresults(j).g.begin_date,1),'.png']);end
  
end

