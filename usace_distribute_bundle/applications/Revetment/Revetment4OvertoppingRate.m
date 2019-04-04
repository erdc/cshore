clear 

% global figID
addpath ../../mfiles

% % load lab data
% % check page 36, 50 in Ali's report
revetment1 = struct ('id', 'RO20B1', 'd1', 38.9/100, 'Tp',  2.3, 'Hrms', 11.3/100, 'dt', 20.6/100, 'Rc', 9.7/100, 'etabar', 0.25/100, 'SWL', 0, 'qosmeasure', 0.28, 'qosali', 3.22) ; 
revetment2 = struct ('id', 'RO20C1', 'd1', 38.9/100, 'Tp',  3.0, 'Hrms', 7.3/100, 'dt', 20.6/100, 'Rc', 9.7/100,  'etabar', 0.37/100, 'SWL', 0, 'qosmeasure', 0.11, 'qosali', 1.41) ; 
revetment3 = struct ('id', 'RO22B1', 'd1', 40.9/100, 'Tp',  2.3, 'Hrms', 11.6/100, 'dt', 22.6/100, 'Rc', 7.7/100, 'etabar', -0.10/100, 'SWL', 2/100, 'qosmeasure', 3.36, 'qosali', 5.02) ; 
revetment4 = struct ('id', 'RO22C1', 'd1', 40.9/100, 'Tp',  2.9, 'Hrms', 7.5/100, 'dt', 22.6/100, 'Rc', 7.7/100, 'etabar', 0.06/100, 'SWL', 2/100, 'qosmeasure', 1.07, 'qosali', 1.73) ; 
revetment5 = struct ('id', 'RO24B1', 'd1', 42.9/100, 'Tp',  2.3, 'Hrms', 11.9/100, 'dt', 24.6/100, 'Rc', 5.7/100, 'etabar', -0.74/100, 'SWL', 4/100, 'qosmeasure', 12.18, 'qosali', 4.53) ; 
revetment6 = struct ('id', 'RO24C1', 'd1', 42.9/100, 'Tp',  2.9, 'Hrms', 7.8/100, 'dt', 24.6/100, 'Rc', 5.7/100, 'etabar', -0.30/100, 'SWL', 4/100, 'qosmeasure', 4.33, 'qosali', 2.38) ; 

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

% lzhu changed here in.iperm to 1, when iperm=1 and in.isedav=1 read zp
% from the makeinfile file, there are variables in.x_p and in.zb_p for the
% profile of the permeable rubble bed.
% Also, when in.iperm = 1, CSHORE read stone diameter, porosity, stability, etc
in.iperm  = 1;          % 0 = no permeability, 1 = permeable

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
in.iveg   = 0;          % vegitation effect
                        % 0: no vegetation or vegetation represented by increased
                        % 1: veg. density, width, height and root depth are 
                        %    specified as input. The height and root depth 
                        %    vary with the bottom elevation change
                        % 2: veg. with constant density, width and height

if in.iveg==3                        
    in.iFv = 1 ; 
    in.idiss  = 2;          % energy dissipation due to vegitation 
                        % (0: veg is accounted in bottom friction (original)
                        %  1: mendez, 2: chen_zhao, 3. use measured wave spectrum)  
end

% in.veg_Cd is moved into the loop
in.veg_n      = 424;       % vegitation density
in.veg_dia    = 0.008;       % vegitation diam
in.veg_ht     = 0.21;          % vegitation height
in.veg_rod    = 0.3;         % vegitation erosion limit below sand for failure
in.veg_extent = [0 0.33]; % vegitation coverage as fraction of total domain length
                                        % [xveg_start, xveg_end] ./ total domain length
in.gamma  = .8;         % shallow water ratio of wave height to water depth
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

%%% Cd curves from Ranjit
cdkc=@(kc) 70.0*kc.^(-0.86);
cdre=@(re) 2*(1300./re + 0.18);
in.veg_Cd = 1.0 ;        % vegitation drag coeff
in.veg_Cdm = in.veg_Cd ;    

in.angle = 0*dum;    % constant incident wave angle at seaward boundary in

in.JONSWAPgamma = 3.3 ;

if in.iveg==3 && in.idiss==2 
    in.numfreq= 500*ones(size(in.freqmax)) ; 
elseif in.iveg==3 && in.idiss==3
    spectrumdata = load ('Jadhav_Omeg_Se.txt') ;
    [row, col] = size(spectrumdata)  ;
    in.numfreq = col *ones(size(in.freqmax)) ; 
end

for icase = 1 : 6
    eval (['in.Hrms = revetment', num2str(icase), '.Hrms *dum;'] )   ; % Hs=sqrt(2)*Hrms, Hrms=Hs/sqrt(2)
    eval (['in.Tp = revetment', num2str(icase), '.Tp*dum ;']) ;  % constant spectral peak period in seconds
    eval (['d1 = revetment', num2str(icase), '.d1 ;']) ; 
    eval (['dt = revetment', num2str(icase), '.dt ;']) ; 
    eval (['Rc = revetment', num2str(icase), '.Rc ;']) ; 
    eval (['swl = revetment', num2str(icase), '.SWL ;']) ; 
    eval (['etabar = revetment', num2str(icase), '.etabar ;']) ; 

    anguphase = 2*pi./in.Tp ; 
    in.freqmin= 0.1*anguphase ; 
    in.freqmax= 8.0*anguphase ; 

    in.Wsetup = etabar*dum ; % wave setup at seaward boundary in meters \bar{eta}
    in.swlbc = swl*dum;    % water level at seaward boundary in meters zb         
%     in.Wsetup = 0*dum ; % wave setup at seaward boundary in meters \bar{eta}
%     in.swlbc = 0*dum;    % water level at seaward boundary in meters zb         
    
    % Idealized numerical tank
    % ZPINP(J,L)= dimensional vertical coordinate in meters of 
    % porous layer bottom or hard or clay bottom at point (J) with 
    % ZPINP(J) equal to or less than ZBINP(J,L) where ZPINP(1,L)=ZBINP(1,L) imposed
    %     zb up, zp down
    Lx       = 8.36;                  % length of domain
    in.dx     = 0.002;       % constant dx 
    in.x_p   = 0:in.dx:Lx;
    in.x       = 0:in.dx:Lx;
    
    x_p         = [0,                     6.3,      6.3+0.82,                                6.3+0.82+1.24] ;
    zb_p     = [-dt-6.3/34.4,       -dt,       -dt+0.82/34.4,                            Rc-0.03] ; 
    in.zb_p = interp1(x_p, zb_p, in.x_p);  
    
    x_b       = [0,                        6.3,      6.3+0.82+1.24-0.55,             6.3+0.82+1.24] ;
    zb         = [-dt-6.3/34.4,        -dt,         Rc,                                         Rc] ; 
    in.zb     = interp1(x_b, zb, in.x); 
    idtmp    = find( isnan(in.zb) == 1); 
    in.zb(idtmp) = in.zb_p(idtmp) ;

    %% friction factor
    fric_fac_smooth = .15;     % bottom friction factor for impermeable bottom
    fric_fac_rough = .5;     % bottom friction factor for permeable cobbles
%     fric_fac_smooth = .002;     % bottom friction factor for impermeable bottom
%     fric_fac_rough = .01;     % bottom friction factor for permeable cobbles
%     
    in.fw     = fric_fac_smooth*ones(size(in.x)); % cross-shore values of bot fric
    idtmp    = find (abs(in.zb-in.zb_p)>1e-10) ; 
    in.fw(idtmp)     = fric_fac_rough; % cross-shore values of bot fric
    
    % note that hp = zb-zp; vertical thickness of the stone layer
% % % % % 
% % % % % figure; hold on; box on
% % % % % plot (in.x, in.zb, 'b')
% % % % % plot (in.x_p, in.zb_p, 'm')

    [kh]  = dispersion(abs(in.zb(1)), in.Tp(1))  ;
    wvlength = 2*pi/(kh/abs(in.zb(1))) ; 
    hv2h  = in.veg_ht /abs(in.zb(1));

    %%% produce input file %%%
    makeinfile_usace(in) ;

    unix(['./../../../src-repo/updatedveg']) ;  
%     results = load_results_usace;
%     unix('rm -f O*') ; 
    
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
    eval (['revetment', num2str(icase), '.qosmodel = qosmodel (icase) ;']) ;

    eval (['qosmeasure(icase) = revetment', num2str(icase), '.qosmeasure ;']) ; 
    eval (['qosali(icase) = revetment', num2str(icase), '.qosali ;']) ; 
end

[qosmodel', qosmeasure']

figure(99); hold on; box on
plot (qosmeasure, qosmodel, '.r', 'markersize', 10, 'linewidth', 2)
plot (qosmeasure, qosali, 'sb', 'markersize', 10, 'linewidth', 2)
fplot (@(x) x, '-k', 'linewidth', 2, 'HandleVisibility', 'off')
fplot (@(x) 1.2*x, ':k', 'linewidth', 2, 'HandleVisibility', 'off')
fplot (@(x) 0.8*x, ':k', 'linewidth', 2, 'HandleVisibility', 'off')
legend ('our results', 'Ali''s results', 'fontsize', 15)
axis equal
xlim ([0, 35])
ylim ([0, 35])
xlabel ('measured q_{os} (cm^2/s)')
ylabel ('modeled q_{os} (cm^2/s)')
set (gca, 'fontsize', 15)

figure(221); hold on; box on
x3 = in.x ;
y31 = in.zb ;
y32 = zeros(size(y31)) ; 
fill ([x3(1:10:end), fliplr(x3(1:10:end))], [y32(1:10:end), fliplr(y31(1:10:end))], rgb('SkyBlue'), 'HandleVisibility', 'off')

% plot (in.x_p, in.zb_p, '-.k', 'linewidth', 2, 'HandleVisibility', 'off')
x3 = in.x ;
y31 = in.zb ;
y32 = in.zb(1)*ones(size(y31)) ; 
fill ([x3(1:10:end), fliplr(x3(1:10:end))], [y32(1:10:end), fliplr(y31(1:10:end))], rgb('SandyBrown'), 'HandleVisibility', 'off')
xlim ([in.x(1), in.x(end)])

% plot (in.x, in.zb, 'b')
% plot (in.x_p, in.zb_p, 'm')

x3 = in.x ;
y31 = in.zb ;
y32 = in.zb_p ; 
fill ([x3(1:10:end), fliplr(x3(1:10:end))], [y32(1:10:end), fliplr(y31(1:10:end))], rgb('SlateGrey'), 'HandleVisibility', 'off')

% plot (in.x, results.hydro.setup, '-b', 'linewidth', 2)
% plot (in.x, results.hydro.Hrms, '-b', 'linewidth', 2)

set (gca, 'fontsize', 15)
xlabel ('x (m)')
ylabel ('elevation (m)')
xlim([0, 10])
ylim([-0.5, 0.4])


% fid1 = fopen ('R2pmodel_revetment.txt', 'w') ; 
% fprintf (fid1, '%f ', R2pmodel) ; 
% fclose (fid1) ; 
% 
% fid2 = fopen ('Htoemodel_revetment.txt', 'w') ; 
% fprintf (fid2, '%f ', H_toe) ; 
% fclose (fid2) ; 

rmpath ../../mfiles

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
