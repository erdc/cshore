%make 5 infiles, centered around times of survey dates
addpath ../../mfiles/
if ~exist('lidar')|~exist('bathy')
  load ./data/bathy.mat
  load ./data/wave.mat
  load ./data/wl.mat

end
waveind = find(strcmp({wave.name},'AWAC-6m'));
g.y = 940;% transects of alts are at y = 770 and 860
[j1 indy] = min(abs(bathy.y_frf-g.y));
g.x_offset = wave(waveind).x_frf;
bathyinds = [70:71]; 

for j = bathyinds
  g.zb_init = bathy.zb(indy,:,j);
  g.zb_fin = bathy.zb(indy,:,j+1);
  g.x_frf = bathy.x_frf;
  % Model simulation timing params
  g.begin_date = bathy.date(j);
  g.dt_bc = 1/24; % day; the time interval of the boundary condition data
  g.end_date = bathy.date(j+1);
  g.bc_dates = g.begin_date:g.dt_bc:g.end_date;
  in.timebc_wave = (g.bc_dates-g.begin_date)*24*3600;
  in.timebc_surg = in.timebc_wave;
  in.nwave = length(in.timebc_wave); in.nsurg = in.nwave;dum = ones(1,in.nwave);

  % CSHORE execution and physical params
  in.header = {'------------------------------------------------------------'
               'CSHORE applied to FRF morphology change'
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
  in.gamma  = .6;      % shallow water ratio of wave height to water depth
  in.sporo  = 0.4;     % sediment porosity                        
  in.d50    = .3;      % d_50 in mm
  in.wf     = vfall(in.d50,20,0); % fall velocity
  in.sg     = 2.65;    % specific gravity
  in.effb   = 0.002;   % suspension efficiency due to breaking eB     
  in.efff   = 0.004;    % suspension efficiency due to friction ef 
  in.slp    = .5;      % suspended load parameter               
  in.slpot  = .1;      % overtopping suspended load parameter               
  in.tanphi = .630;    % tangent (sediment friction angle)        
  in.blp    = 0.001;   % bedload parameter                        
  in.rwh = .01;        % numerical rununp wire height 
  in.ilab = 0;         % controls the boundary condition timing.
  in.fric_fac = .015;  % bottom friction factor

  in.x = 0:in.dx:g.x_offset;
  in.zb = interp1(g.x_frf,g.zb_init,g.x_offset-in.x);
  in.fw = in.fric_fac*ones(size(in.zb)); % cross-shore values of bot fric
  
  d = wave(waveind).date;
  [d ind] = unique(d);
  Hrms = wave(waveind).Hs(ind)/sqrt(2);
  Tm = wave(waveind).Tm(ind);
  Dm_wrtfrf = wave(waveind).Dm_wrtfrf(ind);
  [d2 Hrms] = gap_indicate(d,Hrms,.5);
  in.Hrms = interp1(d2,Hrms,g.bc_dates);
  in.Tp = interp1(d,Tm,g.bc_dates);
  in.angle = interp1(d,Dm_wrtfrf,g.bc_dates);
  indbad = isnan(in.Hrms);
  in.T (indbad ) = NaN;
  in.angle(indbad) = NaN; 
  in.swlbc = interp1(wl.date,wl.wl,g.bc_dates);
  in.Wsetup = 0*in.swlbc;   % wave setup at seaward boundary in meters
  
  %fill in the NaN's with the AWAC-11m
  % waveind11 = find(strcmp({wave.name},'AWAC-11m'));
  % d11 = wave(waveind11).date;
  % [d11 ind] = unique(d);
  % Hrms11 = wave(waveind11).Hs(ind)/sqrt(2);
  % [d11 Hrms11] = gap_indicate(d11,Hrms11,.5);
  %indgood = ~isnan(Hrms11);
  %d11 = d11(indgood);
  %Hrms11 = Hrms11(indgood);
  %[H1,L1,alpha1] = shoal (Hrms11,11,wave(waveind11).Dm_wrtfrf,wave(waveind11).Tm,6);
  %in.Hrms(indbad) = interp1(d11,H1,g.bc_dates(indbad));
  in.filename = sprintf('%03d',j);
  makeinfile_usace(in);
  if exist([in.filename,'.infile'],'file')==2
    save(['g_in',sprintf('%03d',j),'.mat'],'g','in'); 
  end
  
  
end




return


load ./data/summary_bathy.mat;
load ./data/summary_waves_bc.mat;
load ./data/summary_wl.mat;


bnd_gage_num = 1;
for j = 2:3;% bathy ind
  g.begin_date = bathy(j).date-10/24; % start a bit previous to survey
  g.end_date = bathy(j).date+10/24;
  run_inds = find([lidar.date]>=g.begin_date&[lidar.date]<=g.end_date);
  [junk bestind ]=min(abs([lidar.date]-bathy(j).date));

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
  zb = interp1(lidar(bestind).frf_xi,lidar(bestind).zb_full,x_frf);

  in.x = x;
  in.zb = zb;
  in.fw = in.fric_fac*ones(size(in.zb)); % cross-shore values of bot fric

  makeinfile_usace(in);
  movefile('infile',[num2str(j),'.infile']);
  save(['g',num2str(j),'.mat'],'g'); 

end