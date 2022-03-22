iprint =   1;
load results
load ./data/zb.dat
load ./data/x.dat
load ./data/zb_time.dat
cnt = 0;
figure(1);clf;clear hh hlabs
for i = 1:10:length(zb_time)
  cnt = cnt+1;
  [j1 j2] = min(abs(zb_time(i)-[results.morpho.time]));
  hh(cnt)=plot(x(:,i),zb(:,i));hold all
  plot(results.morpho(j2).x,results.morpho(j2).zb,'--','color',get(hh(cnt),'color'));hold all
  hlabs{cnt}= ['Model after  ',num2str(results.morpho(j2).time/3600),' hrs'];
end
legend(hh,hlabs,'location','northwest')
%axis([0 45 -0 7])
ylabel('$z_{b} [m]$','fontname','times','fontsize',20,'interpreter','latex')
xlabel('x [m]','fontname','times','fontsize',14,'fontangle','italic')
title('Profile Evolution','fontname','times','fontsize',14,'fontangle','italic')


if iprint;print -dpng -r300 morph_results.png;end


