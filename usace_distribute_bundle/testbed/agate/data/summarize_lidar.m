data_name = '~/data/agate/data_from_julie_dec_2020/Agate_data_for_Brad_Johnson/';
exp_cases = {'201309292100';'201310161100'}
clear R_all

for i =1:length(exp_cases)
  clear RUNUP
  %if ~exist('RUNUP') 
  %load([data_name,exp_cases{i},'/lidar/20130929-2100.mat'])
  fn = dir([data_name,exp_cases{i},'/lidar/*.mat'])
  load([fn.folder,'/',fn.name])
  %end

  % do my own r2p:
  z = RUNUP.Z_10;
  t = RUNUP.runupT; % I think that this is correct
                    %dt = 24*3600*(t(11)-t(1))/10;
  z_mean = nanmean(z);
  z_nomean = z-z_mean;
  z_nomean_sm = window(z_nomean,7);
  [R] = upcrossing_runup(z_nomean_sm);% does not include mean yet
  R.peaks = R.max_value + z_mean;
  [N,edges]=histcounts(R.peaks,20);
  centers = .5*(edges(2:end)+edges(1:end-1));
  dR = edges(2)-edges(1);
  N_norm = N/(sum(N)*dR);
  cdf = 1-cumsum(N_norm*dR);
  thresh = .02;% for 2% runup
  crsspt = find((cdf(2:end)-thresh).*(cdf(1:end-1)-thresh)<=0);
  crsspt = crsspt(1);
  R.r2p = interp1(cdf(crsspt:crsspt+1),centers(crsspt:crsspt+1),thresh);
  R.r13 = interp_brad(centers,cdf-.333);
  R.mean_runup = z_mean;
  R.std_runup = std(z);
  R.cdf_peaks = cdf;
  R.centers_peaks = centers;
  R.runupT = RUNUP.runupT;
  R.X_10   = RUNUP.X_10  ;
  R.X_20   = RUNUP.X_20  ;
  R.Z_10   = RUNUP.Z_10  ;
  R.Z_20   = RUNUP.Z_20  ;
  R.Z2     = RUNUP.Z2;
  R.x      = RUNUP.x;

  R_all(i).r2p = R.r2p;
  R_all(i).r13 = R.r13;
  R_all(i).mean_runup = R.mean_runup;
  R_all(i).exp_case = exp_case;
end

save R_all R_all 

fs = 26;
figure(1);clf 
%plot(RUNUP.runupT,RUNUP.Z_10,'k','linewidth',2);hold all
plot(t,z_nomean_sm+z_mean,'k','linewidth',2);hold all
datetick
title(['Vertical Runup ',exp_case],'interpreter','latex','fontsize',fs)
xlabel('$time$','interpreter','latex','fontsize',fs)
ylabel('$\eta_r [m]$','interpreter','latex','fontsize',fs)
%print -dpng TD_runup_ts.png
plot(t(R.max_ind),R.peaks,'ro','markerfacecolor','r')
plot(t,R.r2p*ones(size(t)),'r')
%print -dpng TD_runup_ts2.png

ri = R.mean_runup:.01:R.mean_runup+1.5;
cdf2 = exp(-2*((ri-R.mean_runup)/(1.*R.r13-R.mean_runup)).^2); 
figure(2);clf
clear hh
hh(1) = plot(centers,cdf,'b','linewidth',2);hold on
hh(2) = plot(ri,cdf2,'k','linewidth',2);hold on
plot(R.r2p,.02,'ro','markerfacecolor','r')
title(['Runup Peak ',exp_case],'interpreter','latex','fontsize',fs)
xlabel('$\eta_{peak}$','interpreter','latex','fontsize',fs)
ylabel('$cdf$','interpreter','latex','fontsize',fs)
legend(hh,'Data','Rayleigh')
%print -dpng TD_runup_cdf.png