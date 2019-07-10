?clear 

% global figID
addpath ../../mfiles

vegtest = 2 ; % 1 for Test 1 in monthly report, 2 for Test 2 in monthly report

% load lab data
revetment1 = struct ('id', 'RS20B1', 'd1', 38.9/100, 'Tp',  2.3, 'Hrms', 11.2/100, 'dt', 20.6/100, 'Rc', 18.6/100, 'etabar', 0.28/100, 'SWL', 0, 'R2pmeasure', 13.31/100) ; 
revetment2 = struct ('id', 'RS20C1', 'd1', 38.9/100, 'Tp',  3.0, 'Hrms', 7.3/100, 'dt', 20.6/100, 'Rc', 18.6/100, 'etabar', 0.17/100, 'SWL', 0, 'R2pmeasure', 13.59/100) ; 
revetment3 = struct ('id', 'RS22B1', 'd1', 40.9/100, 'Tp',  2.3, 'Hrms', 11.6/100, 'dt', 22.6/100, 'Rc', 16.6/100, 'etabar', 0.11/100, 'SWL', 2/100,  'R2pmeasure', 14.07/100) ; 
revetment4 = struct ('id', 'RS22C1', 'd1', 40.9/100, 'Tp',  2.9, 'Hrms', 7.6/100, 'dt', 22.6/100, 'Rc', 16.6/100, 'etabar', 0.13/100, 'SWL', 2/100,  'R2pmeasure', 14.60/100) ; 
revetment5 = struct ('id', 'RS24B1', 'd1', 42.9/100, 'Tp',  2.3, 'Hrms', 11.9/100, 'dt', 24.6/100, 'Rc', 14.6/100, 'etabar', -0.27/100, 'SWL', 4/100,  'R2pmeasure', 14.12/100) ; 
revetment6 = struct ('id', 'RS24C1', 'd1', 42.9/100, 'Tp',  2.9, 'Hrms', 7.8/100, 'dt', 24.6/100, 'Rc', 14.6/100, 'etabar', 0, 'SWL', 4/100, 'R2pmeasure', 14.98/100) ;   

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
in.iveg   = 3;          % vegitation effect
                        % 0: no vegetation or vegetation represented by increased
                        % 1: veg. density, width, height and root depth are 
                        %    specified as input. The height and root depth 
                        %    vary with the bottom elevation change
                        % 2: veg. with constant density, width and height

if in.iveg==3                        
    in.iFv = 1 ; 
    in.idiss  = 1;          % energy dissipation due to vegitation 
                        % (0: veg is accounted in bottom friction (original)
                        %  1: mendez, 2: chen_zhao, 3. use measured wave spectrum)  
end

in.vegtype   = 2 ; 
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
in.veg_n      = [400, 400];       % vegitation density
in.veg_dia    = [0.01, 0.01];       % vegitation diam
if vegtest == 1
    in.veg_ht     = [0.14, 0.14] ; 
else 
    in.veg_ht     = [0.2, 0.2] ; 
end
in.veg_rod    = [0.3, 0.3];         % vegitation erosion limit below sand for failure
    
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
    
    % Idealized numerical tank
    % ZPINP(J,L)= dimensional vertical coordinate in meters of 
    % porous layer bottom or hard or clay bottom at point (J) with 
    % ZPINP(J) equal to or less than ZBINP(J,L) where ZPINP(1,L)=ZBINP(1,L) imposed
    %     zb up, zp down
    Lx       = 12.0;                  % length of domain
    in.dx     = 0.005;       % constant dx 
    in.x_p   = 0:in.dx:Lx;
    in.x       = 0:in.dx:Lx;

     x_p         = [0,                     6.3,      6.3+0.82,                                6.3+0.82+1.24] ;
    zb_p     = [-dt-0.1831,        -dt,       -dt+0.02384,                            Rc-0.119] ; 
    in.zb_p = interp1(x_p, zb_p, in.x_p, 'linear', 'extrap');  
    
    in.zb     = in.zb_p; 
    
    if vegtest ==1
        in.veg_extent = [6.3/Lx, (6.3+0.82)/Lx; ...
                                  (6.3+0.82)/Lx, (6.3+0.82+1.24)/Lx] ;
    else 
        in.veg_extent = [1.96/Lx, (6.3+0.82)/Lx; ...
                                  (6.3+0.82)/Lx, (6.3+0.82+1.24)/Lx] ;
    end 

    %% friction factor
    fric_fac_smooth = .002;     % bottom friction factor for impermeable bottom
    fric_fac_rough = .01;     % bottom friction factor for permeable cobbles
    
    in.fw     = fric_fac_smooth*ones(size(in.x)); % cross-shore values of bot fric
    idtmp    = find (abs(in.zb-in.zb_p)>1e-5) ; 
    in.fw(idtmp)     = fric_fac_rough; % cross-shore values of bot fric
    
    % note that hp = zb-zp; vertical thickness of the stone layer
% % % % % 
% % % % % figure; hold on; box on
% % % % % plot (in.x, in.zb, 'b')
% % % % % plot (in.x_p, in.zb_p, 'm')

    [kh]  = dispersion(abs(in.zb(1)), in.Tp(1))  ;
    wvlength = 2*pi/(kh/abs(in.zb(1))) ; 
    hv2h  = in.veg_ht /abs(in.zb(1));

    %%% Cd curves from Ranjit
    Awbase           = in.Hrms * cosh(0)/ 2 / sinh(kh) ; % Uw at base of canopy
    omega             = 2*pi / in.Tp(1) ; 
    Uwbase           = Awbase*omega ;
    Rebase            = in.veg_dia * Uwbase / 1e-6 ;
%     KC_jadhav       = Uwbase * (in.Tp/1.35) / 8e-3 ; 
    Cd_jadhav       = .36 + (2600 ./ Rebase) ; 

    in.veg_Cd = Cd_jadhav ;        % vegitation drag coeff
    in.veg_Cdm = in.veg_Cd ;    

    %%% produce input file %%%
    makeinfile_usace_vegfeature(in) ;

    unix(['./../../../src-repo/updatedveg']) ;  
    
    results = load_results_usace;
    
%     unix('rm -f O*') ; 
    
    eval (['revetment', num2str(icase), '.R2pmodel = results.hydro.runup_2_percent ;']) ;
    R2pmodel(icase) = results.hydro.runup_2_percent ;
    eval (['R2pmeasure(icase) = revetment', num2str(icase), '.R2pmeasure ;']) ; 
    
    idtmp = find (in.x>(6.3+0.82)) ; 
    H_toe(icase) = results.hydro.Hrms(idtmp(1)) ;

end

figure(22); hold on; box on; grid on;
load R2pmodel_revetment.txt
if vegtest ==1 
    plot (R2pmodel_revetment, R2pmodel, 'or', 'linewidth', 2, 'markersize', 10)
else 
    plot (R2pmodel_revetment, R2pmodel, '+b', 'linewidth', 2, 'markersize', 10)
end
fplot (@(x) x, '-k', 'linewidth', 2, 'HandleVisibility', 'off')
% fplot (@(x) 1.2*x, ':k', 'linewidth', 2, 'HandleVisibility', 'off')
% fplot (@(x) 0.8*x, ':k', 'linewidth', 2, 'HandleVisibility', 'off')
% text (0.1, 0.07, '-20%', 'fontsize', 15)
% text (0.05, 0.1, '+20%', 'fontsize', 15)
axis equal
xlim ([0, 0.2])
ylim ([0, 0.2])
xlabel ('R_{2%} in Revetment Test Set (m)')
ylabel ('R_{2%} in Vegetation Test Sets (m)')
legend ('Veg. Test Set 1', 'Veg. Test Set 2') 
set (gca, 'fontsize', 15)
set(gca,'linewidth',2)

%%%
figure(33); hold on; box on; grid on;
load Htoemodel_revetment.txt
if vegtest ==1 
    plot (Htoemodel_revetment, H_toe, 'or', 'linewidth', 2, 'markersize', 10)
else 
    plot (Htoemodel_revetment, H_toe, '+b', 'linewidth', 2, 'markersize', 10)
end
fplot (@(x) x, '-k', 'linewidth', 2, 'HandleVisibility', 'off')
axis equal
xlim ([0, 0.2])
ylim ([0, 0.2])
xlabel ('H_{rms} at toe in Revetment Test Set (m)')
ylabel ('H_{rms} at toe in Vegetation Test Sets (m)')
legend ('Veg. Test Set 1', 'Veg. Test Set 2') 
set (gca, 'fontsize', 15)
set(gca,'linewidth',2)

%%%%
figure(44); hold on; box on
x3 = in.x ;
y31 = in.zb ;
y32 = zeros(size(y31)) ; 
fill ([x3(1:10:end), fliplr(x3(1:10:end))], [y32(1:10:end), fliplr(y31(1:10:end))], rgb('SkyBlue'), 'HandleVisibility', 'off')

% plot veg
vegid1 = find (in.x>in.x(end)*in.veg_extent(1) & in.x<in.x(end)*in.veg_extent(2) ) ;
x3 = in.x(vegid1) ;
y31 = in.zb(vegid1) ;
y32 = in.zb(vegid1)+in.veg_ht(1); 
h = fill ([x3(1:10:end), fliplr(x3(1:10:end))], [y32(1:10:end), fliplr(y31(1:10:end))], rgb('green'), 'HandleVisibility', 'off') ;
set(h,'facealpha',.5)

% plot veg
for ii = 1 : in.vegtype
    vegid1 = find (in.x>in.x(end)*in.veg_extent(ii,1) & in.x<in.x(end)*in.veg_extent(ii,2) ) ;
    x3 = in.x(vegid1) ;
    y31 = in.zb(vegid1) ;
    y32 = in.zb(vegid1)+in.veg_ht(ii); 
    h = fill ([x3(1:end), fliplr(x3(1:end))], [y32(1:end), fliplr(y31(1:end))], rgb('green'), 'HandleVisibility', 'off') ;
    set(h,'facealpha',.5)
end

% plot (in.x, in.zb, '-.k', 'linewidth', 2, 'HandleVisibility', 'off')
x3 = in.x ;
y31 = in.zb ;
y32 = in.zb(1)*ones(size(y31)) ; 
fill ([x3(1:10:end), fliplr(x3(1:10:end))], [y32(1:10:end), fliplr(y31(1:10:end))], rgb('SandyBrown'), 'HandleVisibility', 'off')
xlim([0, 10])
set (gca, 'fontsize', 15)
xlabel ('x (m)')
ylabel ('elevation (m)')
xlim([0, 10])
ylim([-0.5, 0.4])

unix('rm -f O*') ; 

rmpath ../../mfiles

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
