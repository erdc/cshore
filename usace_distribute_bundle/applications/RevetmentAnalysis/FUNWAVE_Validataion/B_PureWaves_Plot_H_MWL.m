clear 

fdatadir = '/Users/lzhu/Desktop/OngoingProjects/1_EMI_MomentumBalance/Results/PureWave_T1/scpout/' ;
fdir = '/Users/lzhu/Desktop/OngoingProjects/1_EMI_MomentumBalance/Results/PureWave_T1/' ;

fvcwinclusion = 0 ; 
if fvcwinclusion == 1
    fprintf ('with Fv_cw\n')
else
    fprintf ('no Fv_cw\n')
end

%% Load LWT/NLWT-based H and MWL
load ([fdir, 'PureWaves_H_LWT.txt']) ;
load ([fdir, 'eta_x100m.txt']) ; 
load ([fdir, 'eta_x150m.txt']) ; 

if fvcwinclusion == 0
    PureWaves_MWL_LWT = load ([fdir, 'PureWaves_MWL_LWT_nofvcw.txt']) ;
elseif fvcwinclusion == 1
    PureWaves_MWL_LWT = load ([fdir, 'PureWaves_MWL_LWT_withfvcw.txt']) ;
end

if fvcwinclusion == 0
    PureWaves_MWL_NLWT =  load ([fdir, 'PureWaves_MWL_NLWT_nofvcw.txt']) ;
elseif fvcwinclusion == 1
    PureWaves_MWL_NLWT = load ([fdir, 'PureWaves_MWL_NLWT_withfvcw.txt']) ;
end

H0 = PureWaves_H_LWT(2, 1) ;

m = 2048 ; 
dx = 0.1;
x=[0:m-1]*dx;
xgrid_all = x ; 

% wavemaker and sponge
x_wm = [29.5 30.5 30.5 29.5 29.5];
y_wm = [-5 -5 3 3 -5];
x_sponge_west = [0 20 20 0 0]; 
x_sponge_east = [x(end)-20 x(end) x(end) x(end)-20 x(end)-20]; 
y_sponge_west = [-5 -5 4 4 -5];
y_sponge_east = [-1 -1 4 4 -1];

% bathy
load ([fdatadir, 'depth.txt']) ; 

% for plotting veg and water
vegextend = [130.0, 180.0] ; % in meter 
vegid = find (xgrid_all>vegextend(1) & xgrid_all<vegextend(2)) ; 
x3 = xgrid_all(vegid) ;
y31 = -depth(2, vegid) ;
y32 = -depth(2, vegid) + 4.0 ;
d31 = -depth(2, :) ;
d32 = zeros(size(d31)) ; 

% Plot H
figure; 
subplot (2,1,1); hold on; box on; grid on
h1 = fill ([x, fliplr(x)], [d32(1:end), fliplr(d31(1:end))], rgb('lightblue'), 'HandleVisibility', 'off') ;
h2 = fill ([x3(1:end), fliplr(x3(1:end))], [y32(1:end), fliplr(y31(1:end))], rgb('green'), 'HandleVisibility', 'off') ;
set(h1,'facealpha',.5)
set(h2,'facealpha',.5)
plot (x, -depth, 'k', 'linewidth', 2)
plot(x_wm,y_wm,'b', 'linewidth', 2)
plot(x_sponge_west,y_sponge_west,':k', 'linewidth', 2)
plot(x_sponge_east,y_sponge_east,':k', 'linewidth', 2)
text (3, 3, ['sponge', char(10), 'zone'], 'fontsize', 15)
text (188, 3, ['sponge', char(10), 'zone'], 'fontsize', 15)
text (145, 3.5, 'vegetation zone', 'fontsize', 15)
text (31, 3, ['wave', char(10), 'maker'], 'fontsize', 15)
xlabel('x (m)')
ylabel('bathymetry (m)')
set (gca, 'fontsize', 15)
xlim([0, x(end)])

subplot (2,1, 2); hold on; box on; grid on
for ii = 4
    Havg =load([fdatadir, 'Havg_', num2str(ii, '%05.f')]);
    plot (x, Havg(2, :)/H0,'LineWidth',2)
end
plot (PureWaves_H_LWT(1,:), PureWaves_H_LWT(2, :)/H0,'LineWidth',2)
plot(x_sponge_west, [0, 0, 1.2, 1.2, 0], ':k', 'linewidth', 2, 'HandleVisibility', 'off')
plot(x_sponge_east, [0, 0, 1.2, 1.2, 0], ':k', 'linewidth', 2, 'HandleVisibility', 'off')
legend ('FUNWAVE-TVD', 'theoretical solution based on LWT', 'Location', 'South')
xlabel('x (m)')
ylabel('H/H_0')
set (gca, 'fontsize', 15)
xlim ([0, x(end)])
ylim ([0, 1.2])
set (gcf, 'position', [1849          18        1109         592])

% Plot MWL
load /Users/lzhu/Desktop/OngoingProjects/1_EMI_MomentumBalance/Results/PureWave_T1/MWL.txt
figure; 
subplot (2,1,1); hold on; box on; grid on
h1 = fill ([x, fliplr(x)], [d32(1:end), fliplr(d31(1:end))], rgb('lightblue'), 'HandleVisibility', 'off') ;
h2 = fill ([x3(1:end), fliplr(x3(1:end))], [y32(1:end), fliplr(y31(1:end))], rgb('green'), 'HandleVisibility', 'off') ;
set(h1,'facealpha',.5)
set(h2,'facealpha',.5)
plot (x, -depth, 'k', 'linewidth', 2)
plot(x_wm,y_wm,'b', 'linewidth', 2)
plot(x_sponge_west,y_sponge_west,':k', 'linewidth', 2)
plot(x_sponge_east,y_sponge_east,':k', 'linewidth', 2)
%wave gauges
xgauge1 = [99.5, 100.5] ;
xgauge2 = [149.5, 150.5] ; 
ygauge31 = [1, 1] ;
ygauge32 = [2, 2] ;
h3 = fill ([xgauge1, fliplr(xgauge1)], [ygauge32, fliplr(ygauge31)], rgb('DarkOrange'), 'HandleVisibility', 'off') ;
h4 = fill ([xgauge2, fliplr(xgauge2)], [ygauge32, fliplr(ygauge31)], rgb('DarkOrange'), 'HandleVisibility', 'off') ;
plot ([100, 100], [-0.5, 1], 'color', rgb('DarkOrange'), 'linewidth', 2)
plot ([150, 150], [-0.5, 1], 'color', rgb('DarkOrange'), 'linewidth', 2)
text (3, 3, ['sponge', char(10), 'zone'], 'fontsize', 15)
text (188, 3, ['sponge', char(10), 'zone'], 'fontsize', 15)
text (145, 3.5, 'vegetation zone', 'fontsize', 15)
text (31, 3, ['wave', char(10), 'maker'], 'fontsize', 15)
text (101, 1.5, ['WG1', char(10), 'x=100 m'], 'fontsize', 12)
text (151, 1.5, ['WG2', char(10), 'x=150 m'], 'fontsize', 12)
xlabel('x (m)')
ylabel('bathymetry (m)')
set (gca, 'fontsize', 15)
xlim([0, x(end)])

subplot (2,1, 2); hold on; box on; grid on
plot (MWL (1, :), MWL (2, :)/H0,'LineWidth',2)
plot (PureWaves_MWL_LWT(1,:), PureWaves_MWL_LWT(2, :)/H0, 'LineWidth',2)
plot (PureWaves_MWL_NLWT(1,:), PureWaves_MWL_NLWT(2, :)/H0, 'LineWidth',2)
legend ('FUNWAVE-TVD', 'theoretical solution using LWT-based F_v', 'theoretical solution using SFWT-based F_v', 'location', 'northwest')
xlabel('x (m)')
ylabel('MWL/H_0')
set (gca, 'fontsize', 15)
xlim([0, x(end)])
set (gcf, 'position', [1849          18        1109         592])

% Plot time series of eta at x = 100m and 150m
figure; 
subplot(2,1,1); hold on; box on;
plot (eta_x100m(1, :), eta_x100m(2, :), 'k', 'linewidth', 2)
xlim([49, 201])
ylabel ('\eta at x = 100 m (m)')
set (gca, 'fontsize', 15)

subplot(2,1,2); hold on; box on;
plot (eta_x150m(1, :), eta_x150m(2, :), 'k', 'linewidth', 2)
xlim([49, 201])
xlabel ('time (s)')
ylabel ('\eta at x = 150 m (m)')
set (gca, 'fontsize', 15)
% set (gcf, 'position', [1849          18        1109         592])
