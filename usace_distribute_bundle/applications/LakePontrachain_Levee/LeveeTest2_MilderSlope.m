clear 
addpath ../../mfiles

% This code use small wave heights for case without vegetation. 

% load levee test information
leveeInfo1 = struct ('id', 'case1', 'd1', 9*0.3048, 'Tp',  7.7, 'Hrms', 3.6*0.3048, 'etabar', 0.0, 'SWL', 0.0) ; 

% d1 is the offshore water depth 
figid = 119 ; 
% load Levee profile digitization
LeveeProfile = load ('LeveeProfile_MilderSlope.txt') ;

% prepare for the bathymetry
x1 = LeveeProfile(:, 1)' ;
x_all   = [0, x1] ;
y1 = LeveeProfile(:, 2)' ;
y_all = [LeveeProfile(1,2), y1]  * 0.3048;     % the digitized leveel eleveation is in ft. 1ft = 0.3048 m;
% adjust y_all to d1
y_all = y_all - (y_all(1) + leveeInfo1.d1) ; 


% script params
iplotbc     = 0;                    % 1 to plot the applied boundary conditions
iplotmorpho = 1;                % 1 to plot the computed morphology results
iplothydro  = 1;                  % 1 to plot the computed hydrodynamic results
iplotrunup  = 0;                  % 1 to plot the computed runup position
isave       = 1;                    % 1 to save the computed results

unix('rm -f infile') ;
unix('rm -f O*') ; 

% CSHORE execution and physical params
in.header = {'------------------------------------------------------------'
             'CSHORE applied to idealized planar slope'
             '------------------------------------------------------------'};
in.iline  = 1;          % 1 = single line
in.iprofl = 0;          % 0 = no morph, 1 = run morph
in.isedav = 1;          % 0 = unlimited sand, 1 = hard bottom

in.iperm  = 0;          % 0 = no permeability, 1 = permeable

in.iover  = 1;          % 0 = no overtopping , 1 = include overtopping
in.infilt = 0;          % 1 = include infiltration landward of dune crest
in.iwtran = 0;          % 0 = no standing water landward of crest,
                        % 1 = wave transmission due to overtopping
in.ipond  = 0;          % 0 = no ponding seaward of SWL
in.iwcint = 0;          % 0 = no W & C interaction , 1 = include W & C interaction
in.iroll  = 1;          % 0 = no roller, 1 = roller
in.iwind  = 0;          % 0 = no wind effect
in.itide  = 0;          % 0 = no tidal effect on currents
in.iweibull = 0 ;
in.iveg   = input('in.iveg=') ;          % vegitation effect
                        % 0: no vegetation or vegetation represented by increased
                        % 1: veg. density, width, height and root depth are 
                        %    specified as input. The height and root depth 
                        %    vary with the bottom elevation change
                        % 2: veg. with constant density, width and height
fprintf ('\n') ; 
if in.iveg==3                        
    in.iFv = 2 ; 
    in.idiss  = 1;          % energy dissipation due to vegitation 
                        % (0: veg is accounted in bottom friction (original)
                        %  1: mendez, 2: chen_zhao, 3. use measured wave spectrum)  
end

in.vegtype   = 1 ;
%%%%%
% for in.vegtype =k (k>1)
% in.veg_n = [n1, n2, ..., nk] ; 
% in.veg_dia = [bv1, bv2, ..., bvk] ; 
% same for in.veg_ht, in.veg_rod.
% in.extent = [start1, end1;
%                    start2, end2; 
%                    start3, end3;
%                        :
%                    startk, endk] ;
% in.veg_Cd = [Cd1, Cd2, ..., Cdk]  
% in.veg_Cdm = [Cdm1, Cdm2, ..., Cdmk]
%%%%%
in.veg_n      = 300;       % vegitation density
in.veg_dia    = 0.01;       % vegitation diam
in.veg_ht     = 0.6 ; 
in.veg_rod    = 0.3;         % vegitation erosion limit below sand for failure

in.gamma  = .85;         % shallow water ratio of wave height to water depth
in.sporo  = 0.4;        % sediment porosity                        
in.stoneporo  = 0.5;  % Stone/gravel porosity in porous layer (SNP can be different from sand porosity=0.4 for ISTSAN=1)
in.d50    = 0.1;        % d_50 in mm
in.stonedia = 0.034 ;  % Nominal stone/gravel diameter (m)
in.criticalstability = 1.1 ; % Critical stability number (0.6 to 1.1) for stone
in.wf     = 0.0 ;            %ling changed, was: vfall(in.d50,20,0); % fall velocity
in.sg     = 2.95 ;           % specific gravity
in.effb   = 0.005 ;      % suspension efficiency due to breaking eB     
in.efff   = 0.01 ;       % suspension efficiency due to friction ef 
in.slp    = .5 ;         % suspended load parameter               
in.slpot  = .1;         % overtopping suspended load parameter               
in.tanphi = .630;       % tangent (sediment friction angle)        
in.blp    = 0.001;      % bedload parameter                        
in.rwh    = .02;           % numerical rununp wire height 
in.ilab   = 1;            % controls the boundary condition timing. Don't change

% NOTE: if use in.ilab=0, error message shows up:
% Data input is stopped because the start time for 
% offshore wave conditions and water level is NOT ZERO

% boundary conditions and timing
ftime = 50;         % [sec] final time, dictates model duration
dt = 60;            % time interval in seconds for wave and water level conditions
if in.ilab==0;
  in.timebc_wave = [dt/2:dt:ftime];
else
  in.timebc_wave = [0:dt:ftime];
end
in.timebc_surg = in.timebc_wave;
in.nwave = length(in.timebc_wave); 
in.nsurg = in.nwave; dum = ones(1,in.nwave);

in.angle = 0*dum;    % constant incident wave angle at seaward boundary in

in.JONSWAPgamma = 3.3 ;

if in.iveg==3 && in.idiss==2 
    in.numfreq= 500*ones(size(in.freqmax)) ; 
elseif in.iveg==3 && in.idiss==3
    spectrumdata = load ('Jadhav_Omeg_Se.txt') ;
    [row, col] = size(spectrumdata)  ;
    in.numfreq = col *ones(size(in.freqmax)) ; 
end

for icase = 1
    unix('rm -f O*') ; 
    unix('rm -f *.DAT') ; 
        
    eval (['in.Hrms = leveeInfo', num2str(icase), '.Hrms *dum;'] )   ; % Hs=sqrt(2)*Hrms, Hrms=Hs/sqrt(2)
    eval (['in.Tp = leveeInfo', num2str(icase), '.Tp*dum ;']) ;  % constant spectral peak period in seconds
    eval (['swl = leveeInfo', num2str(icase), '.SWL ;']) ; 
    eval (['etabar = leveeInfo', num2str(icase), '.etabar ;']) ; 

    anguphase = 2*pi./in.Tp ; 
    in.freqmin= 0.05*anguphase ; 
    in.freqmax= 8.0*anguphase ; 

    in.Wsetup = etabar*dum ; % wave setup at seaward boundary in meters \bar{eta}
    in.swlbc = swl*dum;    % water level at seaward boundary in meters zb         
    
    Lx         = 300; % can change to other values (e.g. 310, 350, etc) 
    Nx        = 4000 ;
    in.x       = linspace (0, Lx, Nx);
    in.dx     = in.x(2) - in.x(1) ;       % constant dx 
    
    in.zb      = interp1 (x_all, y_all, in.x, 'linear', 'extrap') ; 
    in.zb_p  = in.zb ; 

    in.vegtype = 1 ; 
    
    in.veg_extent = [245/Lx, (245+45)/Lx] ;    
%     in.veg_extent = [261/Lx, (245+40)/Lx] ;    
                              
    %% friction factor
    fric_fac_smooth = .002;     % bottom friction factor for impermeable bottom
    in.fw     = fric_fac_smooth*ones(size(in.x)); % cross-shore values of bot fric
    
    %%% Cd curves from Jadhav et al. 2013
    [kh]  = dispersion(abs(in.zb(1)), in.Tp(1))  ;
    Awbase           = in.Hrms * cosh(0)/ 2 / sinh(kh) ; % Uw at base of canopy
    omega             = 2*pi / in.Tp(1) ; 
    Uwbase           = Awbase*omega ;
    Rebase            = in.veg_dia * Uwbase / 1e-6 ;
    Cd_jadhav       = .36 + (2600 ./ Rebase) ; 

    in.veg_Cd = Cd_jadhav ;        % vegitation drag coeff
    in.veg_Cdm = in.veg_Cd ;
   
    %%% produce input file %%%
    makeinfile_usace_vegfeature(in) ;

    unix(['./../../../src-repo/Test_fa/testfa']) ;  

    results = load_results_usace;

    %%% load Qos
    fid = fopen('ODOC');
    tot = textscan(fid,'%s','delimiter','\n');
    tot = tot{:};
    fclose(fid);
    
    tmo =strfind(tot,'Total rate (QOTF+QP)')';
    row_ind = find(~cellfun('isempty',tmo));
    row = tot(row_ind);
    col_ind = find(cell2mat(row(1))=='=');
    row = cell2mat(row);
    qosmodel (icase) =str2num(row(:,col_ind+1:end))*1e4;
       
end

figure(figid); hold on; box on; grid on
subplot(311); hold on; box on; grid on
idtmp = find (in.zb(1:end-1).*in.zb(2:end)<0) ; 
x3 = in.x(1:idtmp(1)) ;
y31 = in.zb(1:idtmp(1)) ;
y32 = zeros(size(y31)) ; 
fill ([x3, fliplr(x3)], [y32, fliplr(y31)], rgb('SkyBlue'), 'HandleVisibility', 'off')

x3 = in.x ;
y31 = in.zb ;
y32 = min(in.zb)*ones(size(y31)) ; 
fill ([x3, fliplr(x3)], [y32, fliplr(y31)], rgb('SandyBrown'), 'HandleVisibility', 'off')
xlim ([in.x(1), in.x(end)])

% plot veg
if in.iveg
    for ii = 1 : in.vegtype
        vegid1 = find (in.x>in.x(end)*in.veg_extent(ii,1) & in.x<in.x(end)*in.veg_extent(ii,2) ) ;
        x3 = in.x(vegid1) ;
        y31 = in.zb(vegid1) ;
        y32 = in.zb(vegid1)+in.veg_ht(ii); 
        h = fill ([x3(1:end), fliplr(x3(1:end))], [y32(1:end), fliplr(y31(1:end))], rgb('green'), 'HandleVisibility', 'off') ;
        set(h,'facealpha',.5)
    end
end

x3 = in.x ;
y31 = in.zb ;
y32 = in.zb_p ; 
fill ([x3, fliplr(x3)], [y32, fliplr(y31)], rgb('SlateGrey'), 'HandleVisibility', 'off')

set (gca, 'fontsize', 15)
xlabel ('x (m)')
ylabel ('elevation (m)')
xlim([100, Lx])
ylim([-3, 4])
plot (in.x, in.zb, 'k', 'HandleVisibility', 'off')
if in.iveg==0
    plot (results.hydro.x, results.hydro.Hrms, ':k',  'linewidth', 2)
elseif in.iveg>0
    plot (results.hydro.x, results.hydro.Hrms, ':r',  'linewidth', 2)
end
legend ('H_{rms} w/o veg', 'H_{rms} w/ veg', 'location', 'northwest')

idtmp = find (results.hydro.x>=260) ;
Hrmstoe = results.hydro.Hrms (idtmp(1)) ; 
if in.iveg==0
    text (180, 3.5, ['no veg: H_{rms, 0} = ', num2str(in.Hrms, 2), ' (m)'], 'fontsize', 15) ; 
    text (180, 2.8, ['Qos = ', num2str(qosmodel, 3), ' (cm^2/s)'], 'fontsize', 15) ; 
    text (180, 2.1, ['R2%= ', num2str(results.hydro.runup_2_percent, 2), ' (m)'], 'fontsize', 15) ; 
    text (180, 1.4, ['H_{rms, toe} = ', num2str(Hrmstoe, 2), ' (m)'], 'fontsize', 15) ; 
elseif in.iveg==1
    text (230, 3.5, ['w/ veg (iveg=1): H_{rms, 0} = ', num2str(in.Hrms, 2), ' (m)'], 'fontsize', 15) ; 
    text (230, 2.8, ['Qos = ', num2str(qosmodel, 3), ' (cm^2/s)'], 'fontsize', 15) ; 
    text (230, 2.1, ['R2%= ', num2str(results.hydro.runup_2_percent, 2), ' (m)'], 'fontsize', 15) ; 
    text (230, 1.4, ['H_{rms, toe} = ', num2str(Hrmstoe, 2), ' (m)'], 'fontsize', 15) ; 
elseif in.iveg==3
    text (130, 3.5, ['w/ veg (iveg=3): H_{rms, 0} = ', num2str(in.Hrms, 2), ' (m)'], 'fontsize', 15) ; 
    text (130, 2.8, ['Qos = ', num2str(qosmodel, 3), ' (cm^2/s)'], 'fontsize', 15) ; 
    text (130, 2.1, ['R2%= ', num2str(results.hydro.runup_2_percent, 2), ' (m)'], 'fontsize', 15) ; 
    text (130, 1.4, ['H_{rms, toe} = ', num2str(Hrmstoe, 2), ' (m)'], 'fontsize', 15) ;     
end    
 
subplot(312); hold on; box on; grid on
if in.iveg==0
    plot (results.hydro.x, results.hydro.setup*1000, 'k',  'linewidth', 2)
elseif in.iveg==1
    plot (results.hydro.x, results.hydro.setup*1000, 'r',  'linewidth', 2)
elseif in.iveg==3
    plot (results.hydro.x, results.hydro.setup*1000, 'b',  'linewidth', 2)
end
set (gca, 'fontsize', 15)
% xlabel ('x (m)')
ylabel ('MWL (mm)')
xlim([100, Lx])
ylim ([0, 500])
legend ('MWL w/o veg', 'MWL w/ veg', 'location', 'northwest')

subplot(313); hold on; box on; grid on
if in.iveg==0
    plot (results.hydro.x, results.hydro.umean, ':k', 'linewidth', 2)
elseif in.iveg==1
    plot (results.hydro.x, results.hydro.umean, ':r', 'linewidth', 2)
elseif in.iveg==3
    plot (results.hydro.x, results.hydro.umean, ':b', 'linewidth', 2)
end
set (gca, 'fontsize', 15)
xlabel ('x (m)')
ylabel ('$\overline{u}$ (m/s)', 'interpreter', 'latex')
xlim([100, Lx])
legend ('u_{mean} w/o veg', 'u_{mean} w/ veg', 'location', 'northwest')

set (gcf, 'position', [1485        -179        1391         984])


rmpath ../../mfiles

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
