fs = 12;
close all
if ~exist('R_all')
  load allresults
  load x_offset
  load data/R_all
end

for i = 1:2
  if i==1;exp_case = '201309292100';else;exp_case = '201310161100';end
  fn = ['./data/',exp_case,'/PUV_at_surface/PUV_at_surface.mat'];
  load(fn);
  Hs = 4*std(Pss);
  Hrms = Hs/sqrt(2); 



  figure;clear hh
  hh(1)=plot(allresults(i).results.hydro.x,allresults(i).results.hydro.Hrms,'b-','linewidth',2) ;hold all
  hh(2)=plot(x_offset-xloc_PUV,Hrms,'rs','markerfacecolor','k','markersize',10);
  ylabel('$H_{rms}[m]$','fontname','times','fontsize',fs,'interpreter','latex')
  xlabel('$x[m]$','fontname','times','fontsize',fs,'interpreter','latex')
  set(gca,'TickLabelInterpreter','latex','fontname','times','fontsize',fs)
  title(['Wave Height ',exp_case],'fontsize',fs,'interpreter','latex')
  hl = legend(hh,'CSHORE','Data');set(hl,'fontsize',fs,'interpreter','latex')
  print('-dpng', ['Hrms_',exp_case,'.png']);

  figure;clear hh
  x = allresults(i).results.morpho.x;
  zb = allresults(i).results.morpho.zb;
  eta = allresults(i).results.hydro.setup;
  [xr2p] = interp_brad(x,zb-R_all(i).r2p);
  csr2p = allresults(i).results.hydro.runup_2_percent;
  [csxr2p] = interp_brad(x,zb-csr2p);
  hh(1)=plot(x,eta,'b-','linewidth',2);hold all
  fill([x; x(end); x(1)],[zb; -20; -20],[.9 .9 .7]);hold all
  hh(2)=plot(x_offset-xloc_PUV,mean(Pss)+zloc_PUV_fixed,'rs','markerfacecolor','k','markersize',10);
  hh(3)=plot(csxr2p,csr2p,'bs','markersize',10,'markerfacecolor','b') ;
  hh(4) = plot(xr2p,R_all(i).r2p,'rs','markerfacecolor','r');
  hl = legend(hh,'$\overline{\eta}$ CSHORE','$\overline{\eta}$ Data','$R_{2\%}$ CSHORE','$R_{2\%}$ Data');
  set(hl,'fontsize',fs,'interpreter','latex','location','northwest')
  title(['Mean Free Surface',exp_case],'fontsize',fs,'interpreter','latex')
  ylabel('$\eta [m]$','fontname','times','fontsize',fs,'interpreter','latex')
  xlabel('$x[m]$','fontname','times','fontsize',fs,'interpreter','latex')
  axis([0 1200 -1 6])
  set(gca,'TickLabelInterpreter','latex')
  set(gcf,'Position',[672 545 900 400])
  print('-dpng', ['eta_',exp_case,'.png']);




end





