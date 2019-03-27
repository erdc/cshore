clear 

addpath /Users/lzhu/Desktop/OngoingProjects/5_ERDC/Codes/cshore/src-repo/SurgeReduction_Veg

unix ('rm STRESSES.TXT') ; 
unix ('rm O*') ; 

%% Prepare the bathymetry
load HaihongBahthy.txt
xloc = (70-HaihongBahthy(:, 1))*1000 ;
xloc = xloc - xloc(1) ; 
ydepth = -HaihongBahthy(:,2)  ;
in.dx = 0.4 ; 
in.x = 0:in.dx:xloc(end) ; 
in.zb = interp1 (xloc, ydepth, in.x, 'linear', 'extrap') ; 
initial_zb = in.zb ; 

load HaihongMSL.txt
maxMSL = max(HaihongMSL(:, 2)) ; 
in.zb = in.zb - maxMSL ; 

load HaihongWind.txt
maxwind = max(HaihongWind(:, 2)) ; 

fprintf ('\n\nREMEMBER to change Line531: SURGECOEFF in CSHORE_USACE_VegFeature.f.\n\n')

%%
% Note that manning's n is not the same as bottom friction factor.

% script params
iplotbc     = 0;                  % 1 to plot the applied boundary conditions
iplotmorpho = 1;                  % 1 to plot the computed morphology results
iplothydro  = 1;                  % 1 to plot the computed hydrodynamic results
iplotrunup  = 0;                  % 1 to plot the computed runup position
isave       = 0;                  % 1 to save the computed results

% CSHORE execution and physical params
in.header = {'------------------------------------------------------------'
             'CSHORE applied to idealized planar slope'
             '------------------------------------------------------------'};
in.iline  = 1;          % 1 = single line
in.iprofl = 0;          % 0 = no morph, 1 = run morph
in.isedav = 0;          % 0 = unlimited sand, 1 = hard bottom
in.iperm  = 0;          % 0 = no permeability, 1 = permeable
in.iover  = 1;          % 0 = no overtopping , 1 = include overtopping
in.infilt = 0;          % 1 = include infiltration landward of dune crest
in.iwtran = 0;          % 0 = no standing water landward of crest,
                        % 1 = wave transmission due to overtopping
in.ipond  = 0;          % 0 = no ponding seaward of SWL
in.iwcint = 0;          % 0 = no W & C interaction , 1 = include W & C interaction
in.iroll  = 0;          % 0 = no roller, 1 = roller
in.iwind  = 1;          % 0 = no wind effect
in.itide  = 0;          % 0 = no tidal effect on currents
in.iweibull = 0 ;
in.iveg  = 3 ;          % vegitation effect
in.veg_Cd = 1.0;          % vegitation drag coeff
in.veg_type = 2 ; 

if in.iveg==3                        
    in.idiss  = 1;          % energy dissipation due to vegitation 
                            % (0: veg is accounted in bottom friction (original)
                            %  1: mendez, 2: chen_zhao)  
    in.iFv = 2 ;
    if (in.iFv == 2)
        in.veg_Cdm = in.veg_Cd ;          % vegitation drag coeff
    end
end

% in.veg_extent = [0.4526 0.5458]; % vegitation coverage as fraction of total domain length
in.veg_extent1 = [0.4696, 0.5423]; % vegitation coverage as fraction of total domain length
in.veg_n1  = 255;        % vegitation density
in.veg_dia1 = .0080;        % vegitation diam
in.veg_ht1 = 0.31;        % vegitation height
in.veg_rod1 = .1;         % vegitation erosion limit below sand for failure

in.veg_extent2 = [0.5423, 0.999]; % vegitation coverage as fraction of total domain length
in.veg_n2  = 255;        % vegitation density
in.veg_dia2= .0031;        % vegitation diam
in.veg_ht2 = 1.07;        % vegitation height
in.veg_rod2= .1;         % vegitation erosion limit below sand for failure

in.gamma  = .8;         % shallow water ratio of wave height to water depth
in.sporo  = 0.4;        % sediment porosity                        
in.d50 = .3;            % d_50 in mm
in.wf = 0.0 ; % vfall(in.d50,20,0); % fall velocity
in.sg = 2.65;           % specific gravity
in.effb   = 0.005;      % suspension efficiency due to breaking eB     
in.efff   = 0.01;       % suspension efficiency due to friction ef 
in.slp    = .5;         % suspended load parameter               
in.slpot  = .1;         % overtopping suspended load parameter               
in.tanphi = .630;       % tangent (sediment friction angle)        
in.blp    = 0.001;      % bedload parameter                        
in.rwh = .02;           % numerical rununp wire height 
in.ilab = 1;            % controls the boundary condition timing. Don't change
in.fric_fac = .01;     % bottom friction factor

% boundary conditions and timing
ftime = 1*600;      % [sec] final time, dictates model duration
dt = 1*600;         % time interval in seconds for wave and water level conditions
if in.ilab==1;
  in.timebc_wave = [dt/2:dt:ftime];
else
  in.timebc_wave = [0:dt:ftime];
end
in.timebc_surg = in.timebc_wave;
in.nwave = length(in.timebc_wave); 
in.nsurg = in.nwave;
dum = ones(1,in.nwave);
in.Tp= 2*dum;        % constant spectral peak period in seconds
in.Hrms = 0.001*dum ; 
in.Wsetup = 0*dum;   % wave setup at seaward boundary in meters
in.swlbc = 0*dum; %+1*sin((2*pi/(12*3600))*in.timebc_surg); % water level at seaward boundary in meters
in.angle = 0*dum;    % constant incident wave angle at seaward boundary in

% % Wind parameters
dt_wind = 600 ;
ft_wind  = 1*600 ; 
if in.ilab==1;
    in.time_wind = [dt_wind/2:dt_wind:ft_wind];
else
    in.time_wind = [0:dt_wind:ft_wind];
end

in.nwind = length(in.time_wind) ;
dum = ones(1, in.nwind) ; 
in.wind10  = maxwind * dum ; 
in.windangle = 0*dum ; 
 
% Idealized numerical tank
in.fw = in.fric_fac*ones(size(in.zb)); % cross-shore values of bot fric

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%run cshore%%%%%%%%%%%%%%%%%%%%
makeinfile_usace_varying_vegfeature(in);
unix ('./../cshore_veg') ;

results = load_results_usace;

if results.run_success==1;
  disp(' CSHORE computation completed ');
elseif results.run_success==.5;
  disp(' CSHORE computation partially completed');
end

MWL = results.hydro.setup ; 

%% Plot 
figure(1123); hold on; box on

% Plot bottom bathymetry
plot (in.x/1000, in.zb, '-.k', 'linewidth', 2, 'HandleVisibility', 'off')
% plot (in.x(1:id(1))/1000 , MWL(1:id(1)) , 'color', rgb('darkblue'),
% 'linewidth', 2, 'HandleVisibility','off')
x3 = in.x/1000 ;
y31 = in.zb ;
y32 = in.zb(1)*ones(size(y31)) ; 
fill ([x3(1:50:end), fliplr(x3(1:50:end))], [y32(1:50:end), fliplr(y31(1:50:end))], rgb('SandyBrown'), 'HandleVisibility', 'off')
xlim ([in.x(1)/1000, in.x(end)/1000])

% plot til zero
x3 = in.x/1000 ;
y31 = in.zb ;
y32 = zeros(size(y31)) ; 
fill ([x3(1:50:end), fliplr(x3(1:50:end))], [y32(1:50:end), fliplr(y31(1:50:end))], rgb('SkyBlue'), 'HandleVisibility', 'off')

% plot veg
vegid1 = find (in.x>in.x(end)*in.veg_extent1(1) & in.x<in.x(end)*in.veg_extent1(2) ) ;
x3 = in.x(vegid1)/1000 ;
y31 = in.zb(vegid1) ;
y32 = in.zb(vegid1)+in.veg_ht1; 
h = fill ([x3(1:50:end), fliplr(x3(1:50:end))], [y32(1:50:end), fliplr(y31(1:50:end))], rgb('green'), 'HandleVisibility', 'off') ;
set(h,'facealpha',.5)

vegid2 = find (in.x>in.x(end)*in.veg_extent2(1) & in.x<in.x(end)*in.veg_extent2(2) ) ;
x3 = in.x(vegid2)/1000 ;
y31 = in.zb(vegid2) ;
y32 = in.zb(vegid2)+in.veg_ht2; 
h = fill ([x3(1:50:end), fliplr(x3(1:50:end))], [y32(1:50:end), fliplr(y31(1:50:end))], rgb('green'), 'HandleVisibility', 'off') ;
set(h,'facealpha',.5)

% plot til MWL
if in.iveg==0
    plot (in.x/1000, MWL, '-k', 'linewidth', 2)
elseif in.iveg==1
    plot (in.x/1000, MWL, ':k', 'linewidth', 2)
elseif in.iveg==3
    plot (in.x/1000, MWL, '-.k', 'linewidth', 2)
end
% x3 = in.x/1000 ;
% y31 = zeros(size(y31)) ;
% y32 = MWL' ; 
% fill ([x3(1:50:end), fliplr(x3(1:50:end))], [y32(1:50:end), fliplr(y31(1:50:end))], rgb('SkyBlue'))
% xlim ([in.x(1)/1000, in.x(end)/1000])

xlabel ('x (km)')
ylabel ('elevation & MWL (m)')
set (gca, 'fontsize', 15)
plot (in.veg_extent1(1)*in.x(end)/1000*ones(1,2), [-4, 5], ':k', 'linewidth', 1.5, 'HandleVisibility', 'off')
plot (in.veg_extent1(2)*in.x(end)/1000*ones(1,2), [-4, 5], ':k', 'linewidth', 1.5, 'HandleVisibility', 'off')
plot (in.veg_extent2(2)*in.x(end)/1000*ones(1,2), [-4, 5], ':k', 'linewidth', 1.5, 'HandleVisibility', 'off')
xlim ([in.x(1)/1000, in.x(end)/1000])
ylim([-4, 4.5])
if in.iveg==3
    load MaxMWL_Rigid.txt
    plot ((70-MaxMWL_Rigid(:, 1))-(70-MaxMWL_Rigid(1, 1)) , ...
        MaxMWL_Rigid(:, 2)-MaxMWL_Rigid(1, 2), 'k')
end
set (gca, 'fontsize', 15)
if in.iveg == 3
    legend ('f_0 = 0.035: no veg, surge height = 3.9 m', 'f_0 = 0.035: w/ veg, surge height = -0.5 m', 'Maximum MWL from Zhao and Chen (2016)')
end
% plot original MWL
% plot (in.x/1000, -2.3*ones(size(in.zb)), '--k')

% plot stress
Func_PlotStresses 
load Vel_RigidVeg.txt
plot ((70-Vel_RigidVeg(:, 1))-(70-Vel_RigidVeg(1, 1)) , Vel_RigidVeg(:, 2), 'k')
legend ('U_{surge} (m/s)', 'u_{wave} (m/s)', 'Maximum current speeds (m/s)')

unix ('rm O*') ; 
rmpath /Users/lzhu/Desktop/OngoingProjects/5_ERDC/Codes/cshore/src-repo/SurgeReduction_Veg

