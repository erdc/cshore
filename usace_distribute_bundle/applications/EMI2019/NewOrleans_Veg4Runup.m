clear 

% global figID
addpath ../../mfiles

vegtest = 2 ;

% % load lab data
% % check page 36, 50 in Ali's report
Test1 = struct ('d1', 9*0.3048, 'Tp',  7.7, 'Hrms', 3.6*0.3048, 'Rc', 16.5*0.3048, 'etabar', 0.0, 'SWL', 0.0) ; 

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
in.iveg   = 0;          % vegitation effect
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

in.vegtype = 2; 
in.veg_n      = [400, 400];       % vegitation density
in.veg_dia    = [0.01, 0.01];       % vegitation diam
if vegtest == 1
    in.veg_ht     = [0.14, 0.14] ; 
else 
    in.veg_ht     = [0.2, 0.4] ; 
end
in.veg_rod    = [0.3, 0.3];         % vegitation erosion limit below sand for failure

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
    eval (['in.Hrms = Test', num2str(icase), '.Hrms *dum;'] )   ; % Hs=sqrt(2)*Hrms, Hrms=Hs/sqrt(2)
    eval (['in.Tp = Test', num2str(icase), '.Tp*dum ;']) ;  % constant spectral peak period in seconds
    eval (['d1 = Test', num2str(icase), '.d1 ;']) ; 
    eval (['Rc = Test', num2str(icase), '.Rc ;']) ; 
    eval (['swl = Test', num2str(icase), '.SWL ;']) ; 
    eval (['etabar = Test', num2str(icase), '.etabar ;']) ; 

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
    Lx       = 200+Rc*4;                  % length of domain
    in.dx     = 0.1;       % constant dx 
    in.x_p   = 0:in.dx:Lx;
    in.x       = 0:in.dx:Lx;
    
    x_p         = [0,     200,     200+Rc*4] ;
    zb_p     = [-d1,     -d1,       Rc] ; 
    in.zb_p = interp1(x_p, zb_p, in.x_p);  
    
    in.zb     = in.zb_p ; 
    
    if vegtest ==1
        in.veg_extent = [100/Lx, 150/Lx; ...
                                  150/Lx, 200/Lx] ;
    else 
        in.veg_extent = [100/Lx, 150/Lx; ...
                                  150/Lx, 200/Lx] ;
    end 
    
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
    makeinfile_usace_vegfeature (in) ;

    unix(['./../../../src-repo/updatedveg']) ;  
    results = load_results_usace;
    
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

figure(223); hold on; box on
x3 = in.x ;
y31 = in.zb ;
y32 = zeros(size(y31)) ; 
fill ([x3(1:10:end), fliplr(x3(1:10:end))], [y32(1:10:end), fliplr(y31(1:10:end))], rgb('SkyBlue'), 'HandleVisibility', 'off')

x3 = in.x ;
y31 = in.zb ;
y32 = in.zb(1)*ones(size(y31)) ; 
fill ([x3(1:10:end), fliplr(x3(1:10:end))], [y32(1:10:end), fliplr(y31(1:10:end))], rgb('SandyBrown'), 'HandleVisibility', 'off')
xlim ([in.x(1), in.x(end)])

% plot (in.x, in.zb, 'b')
% plot (in.x_p, in.zb_p, 'm')
% plot veg
for ii = 1 : in.vegtype
    vegid1 = find (in.x>in.x(end)*in.veg_extent(ii,1) & in.x<in.x(end)*in.veg_extent(ii,2) ) ;
    x3 = in.x(vegid1) ;
    y31 = in.zb(vegid1) ;
    y32 = in.zb(vegid1)+in.veg_ht(ii); 
    h = fill ([x3(1:end), fliplr(x3(1:end))], [y32(1:end), fliplr(y31(1:end))], rgb('green'), 'HandleVisibility', 'off') ;
    set(h,'facealpha',.5)
end

x3 = in.x ;
y31 = in.zb ;
y32 = in.zb_p ; 
fill ([x3(1:10:end), fliplr(x3(1:10:end))], [y32(1:10:end), fliplr(y31(1:10:end))], rgb('SlateGrey'), 'HandleVisibility', 'off')

% plot (in.x, results.hydro.setup, '-b', 'linewidth', 2)
% plot (in.x, results.hydro.Hrms, '-b', 'linewidth', 2)

set (gca, 'fontsize', 15)
xlabel ('x (m)')
ylabel ('elevation (m)')
xlim([0, in.x(end)])

rmpath ../../mfiles

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
