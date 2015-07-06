clear all
close all
addpath('../mfiles')

% script params
iplotbc     = 0;                  % 1 to plot the applied boundary conditions
iplotmorpho = 1;                  % 1 to plot the computed morphology results
iplothydro  = 1;                  % 1 to plot the computed hydrodynamic results
iplotrunup  = 1;                  % 1 to plot the computed runup position
isave       = 0;                  % 1 to save the computed results

% CSHORE execution and physical params
in.header = {'------------------------------------------------------------'
             'CSHORE applied to idealized planar slope'
             '------------------------------------------------------------'};
in.iline  = 1;       % 1 = single line
in.iprofl = 1;       % 0 = no morph, 1 = run morph
in.isedav = 0;       % 0 = unlimited sand, 1 = hard bottom
in.iperm  = 0;       % 0 = no permeability, 1 = permeable
in.iover  = 1;       % 0 = no overtopping , 1 = include overtopping
in.infilt = 0;       % 1 = include infiltration landward of dune crest
in.iwtran = 0;       % 0 = no standing water landward of crest,
                     % 1 = wave transmission due to overtopping
in.ipond  = 0;       % 0 = no ponding seaward of SWL
in.iwcint = 0;       % 0 = no W & C interaction , 1 = include W & C interaction
in.iroll  = 0;       % 0 = no roller, 1 = roller
in.iwind  = 0;       % 0 = no wind effect
in.itide  = 0;       % 0 = no tidal effect on currents
in.iveg   = 0;       % vegitation effect
in.dx     = 1;       % constant dx 
in.gamma  = .8;      % shallow water ratio of wave height to water depth
in.sporo  = 0.4;     % sediment porosity                        
in.d50 = .3;         % d_50 in mm
in.wf = vfall(in.d50,20,0); % fall velocity
in.sg = 2.65;        % specific gravity
in.effb   = 0.005;   % suspension efficiency due to breaking eB     
in.efff   = 0.01;    % suspension efficiency due to friction ef 
in.slp    = .5;      % suspended load parameter               
in.slpot  = .1;      % overtopping suspended load parameter               
in.tanphi = .630;    % tangent (sediment friction angle)        
in.blp    = 0.001;   % bedload parameter                        
in.rwh = .02;        % numerical rununp wire height 
in.ilab = 0;         % controls the boundary condition timing. Don't change
in.fric_fac = .015;  % bottom friction factor

% boundary conditions and timing
ftime = 5*3600;      % [sec] final time, dictates model duration
dt = 1*3600;         % time interval in seconds for wave and water level conditions
if in.ilab==1;
  in.timebc_wave = [dt/2:dt:ftime];
else
  in.timebc_wave = [0:dt:ftime];
end
in.timebc_surg = in.timebc_wave;
in.nwave = length(in.timebc_wave); in.nsurg = in.nwave;dum = ones(1,in.nwave);
in.Tp= 8*dum;        % constant spectral peak period in seconds
in.Hrms = 2+.1*cumsum(dum); 
in.Hrms(2:4)=0;
in.Wsetup = 0*dum;   % wave setup at seaward boundary in meters
in.swlbc = 0+1*sin((2*pi/(12*3600))*in.timebc_surg); % water level at seaward boundary in meters
in.angle = 0*dum;    % constant incident wave angle at seaward boundary in

% Idealized numerical tank
Lx = 300;            % length of domain
zb_off = -8;         % offshore bottom position (should be negative)
zb_on = 8;           % onshore bottom position (should be pos)
flat_length = 100;     % length of flat portion at seaward end of
                     % numerical tank
x = [0 flat_length Lx]; % x points
zb = [zb_off zb_off zb_on]; % zb points
in.x = 0:in.dx:Lx;
[j1 j2] = unique(x); 
in.zb = interp1(x(j2),zb(j2),in.x);
in.zb = in.zb + 0 + 0*min(3,exp((in.x-.65*Lx)/20));
in.fw = in.fric_fac*ones(size(in.zb)); % cross-shore values of bot fric

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%run cshore%%%%%%%%%%%%%%%%%%%%
makeinfile_usace(in);

if isunix
  clean
  !../bin/CSHORE_USACE_LINUX.out >scr.txt
else
  !.\bin\cshore_usace_win.out
end


results = load_results_usace;

if results.run_success==1;
  disp(' CSHORE computation completed ');
elseif results.run_success==.5;
  disp(' CSHORE computation partially completed');
end

if isave
  save summary results
end
plot_results_usace

