clear 

% global figID
addpath ../../../mfiles

vegtest = 1 ; % 1 for Test 1 in monthly report, 2 for Test 2 in monthly report

load xgrid_case6.txt
load depth_case6.txt

% load lab data
revetment6 = struct ('id', 'RS24C1', 'd1', 42.9/100, 'Tp',  1/0.4, 'Hrms', 5.2/100, 'dt', 24.6/100, 'Rc', 14.6/100, 'etabar', 0, 'SWL', 0/100, 'R2pmeasure', 14.98/100) ;   

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
in.veg_n      = 400;       % vegitation density
in.veg_dia    = [0.01];       % vegitation diam
if vegtest == 1
    in.veg_ht     = [0.2] ; 
else 
    in.veg_ht     = [0.2] ; 
end
in.veg_rod    = [0.3, 0.3];         % vegitation erosion limit below sand for failure
    
in.gamma  = .9;         % shallow water ratio of wave height to water depth
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

for icase = 6 
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

    Lx       = xgrid_case6(end);                  % length of domain
    in.dx     = 0.0045;       % constant dx 
    in.x       = 0:in.dx:Lx;

    in.zb = interp1(xgrid_case6, -depth_case6(2, :), in.x, 'linear', 'extrap');  
    
    if vegtest ==1
        vegextend = [6.3, (6.3+0.82+1.24)] + 40 ;
    else 
        vegextend = [1.96, (6.3+0.82+1.24)] + 40 ;
    end 
    in.veg_extent  = vegextend / Lx ;

    %% friction factor
    fric_fac_smooth = .002;     % bottom friction factor for impermeable bottom
    fric_fac_rough = .002;     % bottom friction factor for permeable cobbles
    
    in.fw     = fric_fac_smooth*ones(size(in.x)); % cross-shore values of bot fric
    
    % use different fric for vegetation zone
    idtmp    = find ( in.x>in.veg_extent(1)*Lx & in.x<in.veg_extent(end)*Lx ) ; 
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

    in.veg_Cd = 1.8578 ; %Cd_jadhav ;        % vegitation drag coeff
    in.veg_Cdm = in.veg_Cd ;    

    %%% produce input file %%%
    makeinfile_usace_vegfeature(in) ;

    unix(['./../../../../src-repo/updatedveg']) ;  
    
    results = load_results_usace;
    
%     unix('rm -f O*') ; 
    
    eval (['revetment', num2str(icase), '.R2pmodel = results.hydro.runup_2_percent ;']) ;
    R2pmodel(icase) = results.hydro.runup_2_percent ;
    eval (['R2pmeasure(icase) = revetment', num2str(icase), '.R2pmeasure ;']) ; 
    
    idtmp = find (in.x>(6.3+0.82)) ; 
    H_toe(icase) = results.hydro.Hrms(idtmp(1)) ;

end

%%%%
figure; 
subplot(211); hold on; box on
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
xlim([25, 52])
ylim([-0.45, 0.2])
plot (results.hydro.x, results.hydro.Hrms, 'b', 'linewidth', 2)
load Hrms_00008
plot (xgrid_case6, Hrms_00008(2, :), 'k', 'linewidth', 2)
legend ('CSHORE: H_{rms} (m)', 'FUNWAVE-TVD: H_{rms} (m)',  'location', 'southwest')

subplot(212); hold on; box on
plot (results.hydro.x, results.hydro.umean, 'b', 'linewidth', 2)
set (gca, 'fontsize', 15)
xlabel ('x (m)')
ylabel ('$\overline{u}$ (m/s)', 'interpreter', 'latex')
xlim([25, 52])
  
load umean_00008
plot (xgrid_case6, umean_00008(2, :), ':k', 'linewidth', 2)
legend ('CSHORE', 'FUNWAVE-TVD', 'location', 'southwest')


% %%% load Umean
% fid = fopen('OXVELO');
% % tot = textscan(fid,'%f %f %f\n', 'Delimiter', '  ', 'HeaderLines', 1);
% tot = textscan(fid,'%f %f %f\n', 'HeaderLines', 1);
% fclose(fid);
% plot (tot{1}, tot{2}, ':k', 'linewidth', 2)
 
return
figure(929); 
subplot(311); hold on; box on
plot (results.hydro.x, results.hydro.Hrms, 'k', 'linewidth', 2)
set (gca, 'fontsize', 15)
xlabel ('x (m)')
ylabel ('H_{rms} (m)')
xlim([0, Lx])

subplot(313); hold on; box on
plot (results.hydro.x, results.hydro.setup, 'k', 'linewidth', 2)
set (gca, 'fontsize', 15)
xlabel ('x (m)')
ylabel ('MWL (m)')
xlim([25, Lx])

subplot(312); hold on; box on
%%% load Umean
fid = fopen('OXMOME');
% tot = textscan(fid,'%f %f %f\n', 'Delimiter', '  ', 'HeaderLines', 1);
tot = textscan(fid,'%f %f %f\n', 'HeaderLines', 1);
fclose(fid);
Sxx=tot{2} ;
xgrid = tot{1} ; 
%2nd order accurate difference matrices. 
D1 = fdmatrix(xgrid, 1, 2);
dSxxdx = D1 * Sxx ;
plot (xgrid, dSxxdx, 'k', 'linewidth', 2)
set (gca, 'fontsize', 15)
xlabel ('x (m)')
ylabel ('dS_{xx}/dx (N/m^2)')
xlim([25, Lx])


rmpath ../../../mfiles

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
