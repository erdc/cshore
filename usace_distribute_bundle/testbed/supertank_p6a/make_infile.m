addpath ../../mfiles/
load ./data/sb_dat.mat
i = 1;
in.header = {'------------------------------------------------------------'
             ['CSHORE run using ',sbeach_dat(i).title]
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
in.iroll  = 1;       % 0 = no roller, 1 = roller
in.iwind  = 0;       % 0 = no wind effect
in.itide  = 0;       % 0 = no tidal effect on currents
in.iveg   = 0;       % vegitation effect
in.dx     = 1;       % constant dx 
in.gamma  =.7;       % shallow water ratio of wave height to water depth
                     % in.timebc_wave = Hrms_bc(:,1);
                     % in.timebc_surg = in.timebc_wave;
                     % in.nwave  = length(in.timebc_wave); in.nsurg = in.nwave;
                     % in.Hrms   = Hrms_bc(:,2);
                     % in.Tp     = T_bc(:,2);
                     % in.Wsetup = 0*in.timebc_wave;  % wave setup at seaward boundary in meters
                     % in.swlbc  = 0*in.timebc_wave;   % wave setup at seaward boundary in meters
                     % in.angle  = 0*in.timebc_wave;   % incident wave angle at seaward boundary in
                     % in.x      = zb_init(1,1):in.dx:zb_init(end,1); % x points
                     % in.zb     = interp1(zb_init(:,1),zb_init(:,2),in.x); % zb points
cshore_sim_ftime = 60*sbeach_dat(i).ftime;% final time for cshore run in sec
dt_cs = 60*min(sbeach_dat(i).dtelv,sbeach_dat(i).dtwav); % The time step for the boundary conditions
     if dt_cs==0;dt_cs = cshore_sim_ftime/5;end
     in.timebc = dt_cs:dt_cs:cshore_sim_ftime; % BC time series
in.timebc_wave = in.timebc;
in.timebc_surg = in.timebc;
in.nwave = length(in.timebc);in.nsurg = in.nwave;
dum = ones(size(in.timebc));
if sbeach_dat(i).iwave==0
  in.Tp= sbeach_dat(i).T_const*dum;        % spectral peak period in seconds
  in.Hrms = sbeach_dat.H_s_const/sqrt(2)*dum;   % rms waveheightwave;  
elseif sbeach_dat(i).iwave==1
  in.Tp= interp1(sbeach_dat(i).time_wav*60,sbeach_dat(i).T,in.timebc-dt_cs/2,[],'extrap');   % spectral peak period[s]
  in.Hrms = interp1(sbeach_dat(i).time_wav*60,sbeach_dat(i).H_s/sqrt(2),in.timebc-dt_cs/2,[],'extrap');% rms waveheight;  
end
in.Wsetup = 0*dum;  % wave setup at seaward boundary in meters
if sbeach_dat(i).ielv==0
  in.swlbc= sbeach_dat.elv_const*dum;        % 
elseif sbeach_dat(i).ielv==1
  in.swlbc= interp1(sbeach_dat(i).time_elv*60,sbeach_dat(i).elv,in.timebc-dt_cs/2,[],'extrap');      
end
if sbeach_dat(i).iang==0
  in.angle= sbeach_dat.ang_const*dum;        % incident wave angle at seaward boundary 
elseif sbeach_dat(i).iang==1
  in.angle= interp1(sbeach_dat(i).time_ang*60,sbeach_dat(i).ang,in.timebc-dt_cs/2,[],'extrap');      
end

in.x = -sbeach_dat(i).x;in.x = in.x-min(in.x);in.x = flipud(in.x);
in.zb = sbeach_dat(i).z;in.zb = flipud(in.zb);
in.fw     = .01*ones(size(in.zb)); % friction coeff
in.rwh    = .01; 
in.ilab   = 1;
in.sg     = 2.65;
in.d50    = sbeach_dat.d50;
in.wf     = vfall(in.d50,20,0);
in.effb   = .002; % typical = .005 
in.efff   = .003; % typical = .01 
in.slp    = .5;
in.slpot  = .1;
in.tanphi = .66; %typical = .65
in.blp    = .001; %typical .002
makeinfile_usace(in);
