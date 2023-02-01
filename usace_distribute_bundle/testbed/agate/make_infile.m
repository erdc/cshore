addpath ../../mfiles/
load ./data/bathy

for i = 1:2
  if i==1;exp_case = '201309292100';else;exp_case = '201310161100';end
  fn = ['./data/',exp_case,'/PUV_at_surface/PUV_at_surface.mat'];
  load(fn);
  begin_date = datenum(exp_case,'yyyymmddHHMM');
  gaugenum = 8;ind = find(numvec==gaugenum);
  fz=2; dt=1/fz;t=0:dt:(length(Pss(:, 1))-1)*dt;
  eta = Pss(:,ind);eta = eta-mean(eta);
  numsegs = 10;
  stats = find_stats(t, eta, numsegs);

  % CSHORE execution and physical params
  in.filename = exp_case;
  in.header = {'------------------------------------------------------------'
               'CSHORE applied to Agate, OR'
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
  in.iveg   = 0;       % vegitation effect
  in.dx     = 1;       % constant dx 
  in.gamma  = .5;      % shallow water ratio of wave height to water depth
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
  in.rwh = .02;        % numerical rununp wire height 
  in.ilab = 0;         % controls the boundary condition timing.
  in.fric_fac = .015;  % bottom friction factor



  [j1 bathyind] = min(abs(begin_date-[bathy.date]));
  x_offset = xloc_PUV(ind);
  
  in.x  = [0:x_offset-min(bathy(bathyind).xi)];
  in.zb = interp1(x_offset-[bathy(bathyind).xi xloc_PUV(ind)],[bathy(bathyind).zb_sm -10],in.x);
  in.fw = in.fric_fac*ones(size(in.zb)); % cross-shore values of bot fric
  
  in.Hrms = std(eta)*sqrt(8);
  in.Tp = stats.Tp;
  in.angle = 0;
  in.swlbc = mean(Pss(:,gaugenum))+zloc_PUV_fixed(gaugenum);
  in.Wsetup = 0;   % wave setup at seaward boundary in meters
  in.ilab = 1;
  in.nwave = 1;
  in.nsurg = 1;
  in.timebc_wave=0;
  
  makeinfile_usace(in);
  save x_offset x_offset
end

