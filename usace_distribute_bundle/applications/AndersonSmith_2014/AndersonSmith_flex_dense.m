clear

addpath /Users/lzhu/Desktop/OngoingProjects/5_ERDC/Codes/GitCSHORE/Brad/src-repo/CDEPH2CSHORE

cshoreexe = '../flexveg' ;

%vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
rigidity = 'flex' ; % 'flex', 'rigid'
YoungsEstem        = 172.4e6; %80e6 ;
YoungsEblades     = 5e5  ; %5e5 ; %4e8 ; %5e5; %80e6 ;
Cdformula = 'default' ; % 'infile', 'default', 'jadhav', 'moller', 'garzon'
Cdcap = 45.0 ;
nv = [] ; % nv is determined later based on file name.
bv = 6.4e-3 ;
lstem_1pt = 0.415 ;
lblades_1pt = 0.0 ;
nbladeperstem = 0 ; 
bv_blade    = 0.0 ;
tv_thicknessblade = 0.0 ;
segmentnumber = 1 ; %% veg zone (instead of entire domain) is partitioned into # segments.
rhowater       = 1010.5 ; 
gageloc  = [26.9, 27.4, 27.9, 28.5, 29.5, 31.0, 32.7, 34.4, 36.2] ; 
Lx = 63.9;                  % length of domainL

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% define unified formula in case users set in.ivegCd == 0 and want to determine Cd by themselves. 
fCduniversal = @(Re) 0.57+(1546./Re).^1.11 ;

if strcmp(rigidity, 'flex') && (strcmp(Cdformula, 'jadhav') || strcmp(Cdformula, 'moller') || strcmp(Cdformula, 'garzon')) 
    fprintf('The selected formula is only valid for rigid veg. Please change.\n')
    return
end

% initialize some variables
num_Hs = []; obs_Hs = []; 
icount =  1 ;
startinggageid = 1; 
endinggageid = 8; 

% load wave data 
cases=dir(['./USACE_data/*Hrms_MACE.txt']) ;
filename = cases(icount).name 
data     = load (['./USACE_data/', filename]) ;     
filename = strrep(filename, '.txt', '') ; 
id0      = strfind(filename, '-') ;
str1     = filename(1:id0(1)-1) ;
str2     = filename((id0(1)+1):id0(2)-1) ;
id3      = strfind(filename, 'aq_') ;
waterdepth = str2num([str1(1:end-1)]) * 0.0254 ;  % 1 inch = 0.0254 m
waveperiod = str2num(str2(2:end)) ;
if waveperiod > 100
    waveperiod = waveperiod/100 ; 
else
    waveperiod = waveperiod/10 ; 
end

Tmeanfilename = strrep(filename, 'Hrms', 'Tmean') ;
Tmean = load (['./USACE_data/', Tmeanfilename, '.txt']) ;

if sum(strfind([filename], 'med'))~=0
    nv = 200 ;
elseif sum(strfind([filename], 'high'))~=0
    nv = 400 ;
else
    nv = 0.0 ;
end

%%% load measured data %%%
exp_H = data(5:end-2) ; % this is Hrms

%% start CSHORE inputs
% script params
iplotbc     = 0;                  % 1 to plot the applied boundary conditions
iplotmorpho = 1;                  % 1 to plot the computed morphology results
iplothydro  = 1;                  % 1 to plot the computed hydrodynamic results
iplotrunup  = 0;                  % 1 to plot the computed runup position
isave       = 1;                  % 1 to save the computed results

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

in.dx     = 0.02;       % constant dx 
in.gamma  = 0.7; % .8;         % shallow water ratio of wave height to water depth
in.sporo  = 0.4;        % sediment porosity                        
in.d50    = .3;            % d_50 in mm
in.wf = 0.0;            %ling changed, was: vfall(in.d50,20,0); % fall velocity
in.sg = 2.65;           % specific gravity
in.effb   = 0.005;      % suspension efficiency due to breaking eB     
in.efff   = 0.01;       % suspension efficiency due to friction ef 
in.slp    = .5;         % suspended load parameter               
in.slpot  = .1;         % overtopping suspended load parameter               
in.tanphi = .630;       % tangent (sediment friction angle)        
in.blp    = 0.001;      % bedload parameter                        
in.rwh    = .02;           % numerical rununp wire height 
in.ilab   = 1;            % controls the boundary condition timing. Don't change

% NOTE: if use in.ilab=0, error message shows up:
% Data input is stopped because the start time for 
% offshore wave conditions and water level is NOT ZERO

in.fric_fac = .16;     % bottom friction factor

% boundary conditions and timing
ftime = 50;         % [sec] final time, dictates model duration
dt = 60;            % time interval in seconds for wave and water level conditions
if in.ilab==0
  in.timebc_wave = [dt/2:dt:ftime];
else
  in.timebc_wave = [0:dt:ftime];
end
in.timebc_surg = in.timebc_wave;
in.nwave = length(in.timebc_wave); 
in.nsurg = in.nwave; dum = ones(1,in.nwave);

in.Tp= Tmean*dum ; % constant spectral peak period in seconds
in.Hrms = 0.5*(exp_H(startinggageid)+exp_H(startinggageid+1)) * dum;  % wave height at G5 is considered as
                        % incident wave height
                        % Hrms = Hs * 0.707 ; 
                        % we are using Hrms data directly here!

in.Wsetup = 0*dum ; %%labdata{6}(2)*1e-3*dum ;  % wave setup at seaward boundary in meters \bar{eta}
in.swlbc = 0*dum;    % water level at seaward boundary in meters zb         

in.angle = 0*dum;    % constant incident wave angle at seaward boundary in

% Idealized numerical tank
zb_off = -waterdepth;       % offshore bottom position (should be negative)
zb_on = 1.3-waterdepth;                 % onshore bottom position (should be pos)

flat_length = 37.9;         % length of flat portion at seaward end of numerical tank
x = [0 flat_length Lx];     % x points
zb = [zb_off zb_off zb_on]; % zb points
in.x = 0:in.dx:Lx;
[j1, j2] = unique(x); 

in.zb = interp1(x(j2),zb(j2),in.x);
in.fw = in.fric_fac*ones(size(in.zb)); % cross-shore values of bot fric

%%% veg info
in.iveg   = 3;          % vegitation effect
                        % 0: no vegetation or vegetation represented by increased
                        % 1: veg. density, width, height and root depth are 
                        %    specified as input. The height and root depth 
                        %    vary with the bottom elevation change
                        % 2: veg. with constant density, width and height
in.idiss  = 1  ;        % energy dissipation due to vegitation 
                        % (0: veg is accounted in bottom friction (original)
                        %  1: mendez, 2: chen_zhao)  
in.ibreaking = 0 ; % 1: turn on breaking model (a>=1); 0: turn off breaking model (no a>=1 limitation)
                            % default is 1                           
in.iFv = 1; 

if strcmp(rigidity, 'flex')
    in.ivegtype = 1;
else
    in.ivegtype = 0 ;
end

if strcmp(Cdformula, 'infile')
    in.ivegCd = 0;
elseif strcmp(Cdformula, 'default')
    in.ivegCd = 1;
elseif strcmp(Cdformula, 'jadhav')
    in.ivegCd = 2;
elseif strcmp(Cdformula, 'garzon')
    in.ivegCd = 3;
elseif strcmp(Cdformula, 'moller')
    in.ivegCd = 4;
end

in.nvegsegment = segmentnumber;

in.veg_extent = [ gageloc(startinggageid)- gageloc(startinggageid), gageloc(endinggageid)- gageloc(startinggageid)] / Lx; % vegitation coverage as fraction of total domain length
                        % [xveg_start, xveg_end] ./ total domain length

in.veg_n = nv*ones(size(in.x));
in.veg_dia= bv*ones(size(in.x));

in.veg_nblade      = nbladeperstem *in.veg_n;       % vegitation density
in.veg_diablade    = bv_blade*ones(size(in.x)) ;        % vegitation diam
in.veg_thicknessblade = tv_thicknessblade*ones(size(in.x)) ;  
lstem = lstem_1pt * ones(size(in.x)) ; 
lblades = lblades_1pt * ones(size(in.x)) ; 
in.veg_ht     =  lstem ;          % vegitation height
in.veg_htblade  = lblades ;          % vegitation height          
in.veg_Estem = YoungsEstem * ones(size(in.x)) ; 
in.veg_Eblade = YoungsEblades * ones(size(in.x)) ; 
in.veg_rod    = 0.3* ones(size(in.x))  ;         % vegitation erosion limit below sand for failure
in.rhowater       = rhowater ; 

idnoveg_tmp      = find(in.x < in.veg_extent(1)*Lx | in.x > in.veg_extent(end)*Lx) ;    %%  
in.veg_n(idnoveg_tmp) =  0.0; in.veg_dia(idnoveg_tmp) =  0.0; 
in.veg_ht(idnoveg_tmp) =  0.0;in.veg_Estem(idnoveg_tmp) =  0.0;

in.veg_nblade(idnoveg_tmp) =  0.0; in.veg_diablade(idnoveg_tmp) =  0.0; 
in.veg_thicknessblade(idnoveg_tmp) =  0.0; in.veg_htblade(idnoveg_tmp) =  0.0;
in.veg_Eblade(idnoveg_tmp) =  0.0;

if in.ivegCd == 0 
    %% Determine Cd
    omega        = 2.0*pi/ (Tmean) ;
    id1_tmp      = find(in.x >= (in.veg_extent(1)*Lx)) ;    
    id1          = id1_tmp(1) ; 
    [kh]         = dispersionLZ(abs(in.zb(id1)), Tmean)  ;

    wavenum    	 = kh/abs(in.zb(id1)) ;

    Awbase       = in.Hrms * cosh(0)/ 2 / sinh(kh) ; % Uw at base of canopy
    Uwbase       = Awbase*omega ;
    Rebase   = in.veg_dia(1) * Uwbase / 1e-6 ;

    in.veg_Cd(id1_tmp) = fCduniversal(Rebase) ; 
    in.veg_Cdm = in.veg_Cd ; 

end    

%% Run CSHORE
makeinfile_usace_flexveg_CdEPH2CSHORE(in) ;
unix(cshoreexe) ;  
results = load_results_usace_CdEPH2CSHORE;

%% Find out numerical H
num_H = [] ;
for ii  = 1:length(gageloc)
   id1         = find(results.hydro.x >= gageloc(ii)-gageloc(startinggageid)) ;
   num_H(ii)   = results.hydro.Hrms(id1(1)) ; % this is Hrms
end
num_Hs    = [num_Hs, num_H([startinggageid:endinggageid])] ; 
obs_Hs    = [obs_Hs, exp_H([startinggageid:endinggageid])] ;     


figure; hold on; box on 
plot(in.x, results.hydro.Hrms, 'k')
plot(gageloc([startinggageid:endinggageid])-gageloc(startinggageid), exp_H([startinggageid:endinggageid]), 'or')
xlabel('x (m)')
ylabel('H_{rms} (m)')
xlim([0, 9.8])
ylim([0, inf])
set (gca, 'fontsize', 18)
set(gcf,'color','w');


figure; hold on; box on
plot (obs_Hs, num_Hs, '.g', 'markersize', 20)
fplot (@(x) x, 'k', 'linewidth', 2, 'handlevisibility', 'off')
axis square
xlabel ('obs. Hrms (m)')
ylabel ('num. Hrms (m)')
xlim ([0, 0.2])
ylim ([0, 0.2])
set (gca, 'fontsize', 18)
set(gcf,'color','w');

rmpath /Users/lzhu/Desktop/OngoingProjects/5_ERDC/Codes/GitCSHORE/Brad/src-repo/CDEPH2CSHORE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
