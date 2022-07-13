%make 5 infiles, centered around times of survey dates
addpath ../../mfiles/
if ~exist('lidar_sm')|~exist('bathy')
  load ./data/summary_bathy.mat;
  load ./data/summary_waves_bc.mat;
  load ./data/summary_wl.mat;
  load ./data/summary_lidar_science_small.mat
end
bnd_gage_num = 1;
for j = 2:3;% bathy ind
  g.begin_date = bathy(j).date-10/24; % start a bit previous to survey
  g.end_date = bathy(j).date+10/24;
  run_inds = find([lidar_sm.date]>=g.begin_date&[lidar_sm.date]<=g.end_date);
  [junk bestind ]=min(abs([lidar_sm.date]-bathy(j).date));

  % CSHORE execution and physical params
  in.header = {'------------------------------------------------------------'
               'CSHORE applied to FRF cross-shore array for runup'
               '------------------------------------------------------------'};
  in.iline  = 1;       % 1 = single line
  in.iprofl = 0;       % 0 = no morph, 1 = run morph
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
  in.iveg   = 0;          % vegitation effect
  in.dx     = 1;       % constant dx 
  in.gamma  = .8;      % shallow water ratio of wave height to water depth
  in.sporo  = 0.4;     % sediment porosity                        
  in.d50    = .3;         % d_50 in mm
  in.wf     = vfall(in.d50,20,0); % fall velocity
  in.sg     = 2.65;        % specific gravity
  in.effb   = 0.005;   % suspension efficiency due to breaking eB     
  in.efff   = 0.01;    % suspension efficiency due to friction ef 
  in.slp    = .5;      % suspended load parameter               
  in.slpot  = .1;      % overtopping suspended load parameter               
  in.tanphi = .630;    % tangent (sediment friction angle)        
  in.blp    = 0.001;   % bedload parameter                        
  in.rwh = .01;        % numerical rununp wire height 
  in.ilab = 1;         % controls the boundary condition timing.
  in.fric_fac = .015;  % bottom friction factor

  ftime = 24*3600*(g.end_date-g.begin_date); % [sec] final time, dictates model duration
  dt = 3600;
  in.timebc_wave = [0:dt:ftime];
  in.timebc_surg = in.timebc_wave;
  in.nwave = length(in.timebc_wave); in.nsurg = in.nwave;dum = ones(1,in.nwave);
  g.date_bc = g.begin_date+in.timebc_wave/(24*3600);
  g.Tm = interp1(wave(bnd_gage_num).time,wave(bnd_gage_num).Tm,g.date_bc);
  g.Hrms = interp1(wave(bnd_gage_num).time,wave(bnd_gage_num).Hs./sqrt(2),g.date_bc);
  g.angle = interp1(wave(bnd_gage_num).time,wave(bnd_gage_num).Dm_wrtfrf,g.date_bc);
  g.wl = interp1(wl.date,wl.wl,g.date_bc);
  in.Tp= g.Tm;        % constant spectral peak period in seconds
  in.Hrms = g.Hrms;
  in.Wsetup = 0*dum;   % wave setup at seaward boundary in meters
  in.swlbc = g.wl; % water level at seaward boundary in meters
  in.angle = g.angle;    % constant incident wave angle at seaward boundary in

  g.x_offset =  wave(bnd_gage_num).frf_x;
  x_frf = g.x_offset:-in.dx:70;
  x = g.x_offset-x_frf;
  zb = interp1(lidar_sm(bestind).frf_xi,lidar_sm(bestind).zb_full,x_frf);

  in.x = x;
  in.zb = zb;
  in.fw = in.fric_fac*ones(size(in.zb)); % cross-shore values of bot fric

  makeinfile_usace(in);
  movefile('infile',[num2str(j),'.infile']);
  save(['g',num2str(j),'.mat'],'g'); 

end