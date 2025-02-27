addpath ../../mfiles/
load ./data/lstf_data.mat
g.x_offset = 25 ;


% CSHORE execution and physical params
in.header = {'------------------------------------------------------------'
             'CSHORE LSTF Test 1 Case 3 '
             '------------------------------------------------------------'};
in.iline  = 1;          % 1 = single line
in.iprofl = 1;          % 0 = no morph, 1 = run morph, 1.1 = run morph without initial smoothing
in.isedav = 0;          % 0 = unlimited sand, 1 = hard bottom
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
in.iveg   = 0;          % vegitation effect
in.veg_Cd = 1;          % vegitation drag coeff
in.veg_n  = 100;        % vegitation density
in.veg_dia= .01;        % vegitation diam
in.veg_ht = .20;        % vegitation height
in.veg_rod= .1;         % vegitation erosion limit below sand for failure
in.veg_extent = [.7 1]; % vegitation coverage as fraction of total domain length
in.dx     = .5;          % constant dx 
in.gamma  = 1;         % shallow water ratio of wave height to water depth
in.sporo  = 0.4;        % sediment porosity                        
in.d50 = .15;            % d_50 in mm
in.wf = vfall(in.d50,20,0); % fall velocity
in.sg = 2.65;           % specific gravity
in.effb   = 0.002;      % suspension efficiency due to breaking eB     
in.efff   = 0.004;       % suspension efficiency due to friction ef 
in.slp    = .5;         % suspended load parameter               
in.slpot  = .1;         % overtopping suspended load parameter               
in.tanphi = .630;       % tangent (sediment friction angle)        
in.blp    = 0.003;      % bedload parameter                        
in.rwh = .02;           % numerical rununp wire height 
in.ilab = 1;            % controls the boundary condition timing. Don't change
in.fric_fac = .02;     % bottom friction factor

% boundary conditions and timing
ftime = 3*3600;      % [sec] final time, dictates model duration
dt = 1*3600;         % time interval in seconds for wave and water level conditions
if in.ilab==1;
  in.timebc_wave = [dt/2:dt:ftime];
else
  in.timebc_wave = [0:dt:ftime];
end
in.timebc_surg = in.timebc_wave;
in.nwave = length(in.timebc_wave); in.nsurg = in.nwave;dum = ones(1,in.nwave);
in.Tp= 1.5*dum;        % constant spectral peak period in seconds
in.Hrms = .19*dum; 
in.Wsetup = 0*dum;   % wave setup at seaward boundary in meters
in.swlbc = 0*dum; % water level at seaward boundary in meters
in.angle = 10*dum;    % constant incident wave angle at seaward boundary in

x = g.x_offset-[0;lstf.x_bathy;g.x_offset]; % x points
zb = [lstf.zb_bathy(1);lstf.zb_bathy;lstf.zb_bathy(end)]; % zb points
in.x = min(x):in.dx:max(x);
in.zb = interp1(x,zb,in.x);
in.fw = in.fric_fac*ones(size(in.zb)); % cross-shore values of bot fric

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%run cshore%%%%%%%%%%%%%%%%%%%%
makeinfile_usace(in);
save g_in g in


