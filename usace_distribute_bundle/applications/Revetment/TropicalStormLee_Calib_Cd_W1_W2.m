clear 

% global figID
addpath ../../mfiles

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
% lzhu changed here in.iperm, when iperm=1 and in.isedav=1 read zp
% from the makeinfile file, there are variables in.x_p and in.zb_p for the
% profile of the permeable rubble bed.
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
in.dx     = 0.1;       % constant dx 
in.gamma  = .8;         % shallow water ratio of wave height to water depth
in.sporo  = 0.5;        % sediment porosity                        
in.d50    = 34;            % d_50 in mm
in.wf     = 0.0;            %ling changed, was: vfall(in.d50,20,0); % fall velocity
in.sg     = 2.95;           % specific gravity
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

in.fric_fac = .015;     % bottom friction factor

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

in.Hrms = 0.112*dum; % Hs=sqrt(2)*Hrms, Hrms=Hs/sqrt(2)
in.Wsetup = 0*dum ; %%labdata{6}(2)*1e-3*dum ;  % wave setup at seaward boundary in meters \bar{eta}
in.swlbc = 0*dum;    % water level at seaward boundary in meters zb         
in.angle = 0*dum;    % constant incident wave angle at seaward boundary in

in.Tp= 2.3*dum ; % constant spectral peak period in seconds
anguphase = 2*pi./in.Tp ; 
in.freqmin= 0.1*anguphase ; 
in.freqmax= 8.0*anguphase ; 
in.JONSWAPgamma = 3.3 ;

if in.iveg==3 && in.idiss==2 
    in.numfreq= 500*ones(size(in.freqmax)) ; 
elseif in.iveg==3 && in.idiss==3
    spectrumdata = load ('Jadhav_Omeg_Se.txt') ;
    [row, col] = size(spectrumdata)  ;
    in.numfreq = col *ones(size(in.freqmax)) ; 
end

% Idealized numerical tank
Lx   = 33;                  % length of domain
in.x = 0:in.dx:Lx;
x1   = [0,       6.3+0.7,    33] ;
zb1  = [-0.21, -0.18, 0.2] ; 
in.zb = interp1(x1, zb1, in.x); 

in.x_p = 0:in.dx:Lx;
x2       = [0,       6.3,    33] ;
zb2  = [-0.21, -0.15, 0.2] ; 
in.zb_p = interp1(x2, zb2, in.x_p); 

%%%
in.fw = in.fric_fac*ones(size(in.zb)); % cross-shore values of bot fric

[kh]  = dispersion(abs(zb1(1)), in.Tp(1))  ;
wvlength = 2*pi/(kh/abs(zb1(1))) ; 
hv2h  = in.veg_ht /abs(zb1(1));

in.veg_Cd = 1.0 ;        % vegitation drag coeff
in.veg_Cdm = in.veg_Cd ;    

%%% produce input file %%%
makeinfile_usace(in) ;

return
unix(['./../../../src-repo/updatedveg']) ;  
results = load_results_usace;

unix('rm -f O*') ; 
    
rmpath ../../mfiles

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
