iprint = 1;
load ./data/lstf_data.mat
load ./data/summary_sedflux.mat
load results 
load g_in
clear hh hlabs
fs = 14;

figure;
plot(results.hydro(end).x,results.hydro(end).Hrms,'linewidth',2);hold all
title('$RMS$ Wave Height','fontsize',fs,'interpreter','latex')
ylabel('$H_{rms}[m]$','fontname','times','fontsize',fs,'interpreter','latex')
xlabel('$x[m]$','fontname','times','fontsize',fs,'interpreter','latex')
plot(g.x_offset-lstf.x_waves,lstf.Hrms,'rs','markerfacecolor','k')
if iprint;print -dpng wave_height.png;end  
figure
plot(results.hydro(end).x,results.hydro(end).vmean,'linewidth',2);hold all
title('Longshore Velocity','fontsize',fs,'interpreter','latex')
ylabel('$V[m/s]$','fontname','times','fontsize',fs,'interpreter','latex')
xlabel('$x[m]$','fontname','times','fontsize',fs,'interpreter','latex')
plot(g.x_offset-lstf.x_currents,-lstf.V/100,'rs','markerfacecolor','k')
if iprint;print -dpng longhshore_current.png;end


figure
x_qy = g.x_offset-sedflux.trap_x;
qy = sedflux.flux_corrected_m2_per_s;
plot(results.sed(end).x_cross,results.sed(end).qy,'b','linewidth',2);hold all
plot(x_qy,qy,'rs','markerfacecolor','k')
title('Longshore Sediment Transport ',...
      'fontname','times','fontsize',fs,'interpreter','latex')
ylabel('$Q_y[m^2/s]$','fontname','times','fontsize',fs,'interpreter','latex')
xlabel('$x[m]$','fontname','times','fontsize',fs,'interpreter','latex')
if iprint;print -dpng longhshore_tranport.png;end

