fs = 12;
iprint =   1;
load results
load ./data/sb_dat.mat

figure(1);clf;clear hh hlabs

hh(1)=plot(-(sbeach_dat.x-sbeach_dat.x(end)),sbeach_dat.z,'k','linewidth',2);hold all
hh(2)=plot(-(sbeach_dat.x_final-sbeach_dat.x(end)),sbeach_dat.z_final,'r-','linewidth',2);
hlabs{1}=['Measured initial'];
hlabs{2}=['Measured, t = ',num2str(60*sbeach_dat.ftime),' s'];
hh(3)=plot(results.morpho(end).x,results.morpho(end).zb,'b--','linewidth',2);
hlabs{3}=['CSHORE, t = ',num2str(results.morpho(end).time),' s'];
legend(hh,hlabs,'location','northwest')
title([sbeach_dat.name,' Profile Evolution '],'fontname','times','fontsize',fs,'fontangle','italic')
xlabel('$x[m]$','fontname','times','fontsize',fs,'interpreter','latex')
ylabel('$z_b[m]$','fontname','times','fontsize',fs,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex','fontname','times','fontsize',fs)


if iprint;print -dpng -r300 morph_results.png;end


