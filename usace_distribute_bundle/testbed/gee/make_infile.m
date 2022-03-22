addpath ../../mfiles/

load ./data/Hrms_bc.dat
load ./data/T_bc.dat
load ./data/zb_init.dat

in.header = {'------------------------------------------------------------'
            ['CSHORE Great Erosion Experiment                             ']
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
in.timebc_wave = Hrms_bc(:,1);
in.timebc_surg = in.timebc_wave;
in.nwave  = length(in.timebc_wave); in.nsurg = in.nwave;
in.Hrms   = Hrms_bc(:,2);
in.Tp     = T_bc(:,2);
in.Wsetup = 0*in.timebc_wave;  % wave setup at seaward boundary in meters
in.swlbc  = 0*in.timebc_wave;   % wave setup at seaward boundary in meters
in.angle  = 0*in.timebc_wave;   % incident wave angle at seaward boundary in
in.x      = zb_init(1,1):in.dx:zb_init(end,1); % x points
in.zb     = interp1(zb_init(:,1),zb_init(:,2),in.x); % zb points
in.fw     = .01*ones(size(in.zb)); % friction coeff
in.rwh    = .01; 
in.ilab   = 0;
in.wf     = .03;
in.sg     = 2.65;
in.d50    = .22;
in.effb   = .003; % typical = .005 
in.efff   = .006; % typical = .01 
in.slp    = .5;
in.slpot  = .1;
in.tanphi = .66; %typical = .65
in.blp    = .001; %typical .002

makeinfile_usace(in);

save g_in in




return


data_dir = '/home/johnson/data/crosstex/gee/summary/';
boundary_gauge=1
if ~exist('waves_short')
  load ([data_dir 'waves_short'])  
  load ([data_dir 'profile'])  
end
Hrms_data = [ 0.28 0.33 0.33  0.39  0.42];  
T_data = [4 4 6 4 6];
repeats = [6 18 18 15 64+23];
Hrms = [];
T = [];
for i = 1:length(Hrms_data)
  Hrms = [Hrms ;repmat(Hrms_data(i),repeats(i),1)]
  T = [T ;repmat(T_data(i),repeats(i),1)]
end

offset = ft2m(10+3);
in.header = {'------------------------------------------------------------'
            ['CSHORE of GEE,                                    July, 2007']
             '------------------------------------------------------------'};
in.iprofl = 1;  % 0 = no morph, 1 = run morph
in.iroll = 1;   % 0 = no roller, 1 = roller
in.iwind = 0;   % 0 = no wind effect
in.gamma = .5;  % shallow water ratio of wave height to water depth
dt = 1/4;
in.nwave=36/dt
dum = ones(1,in.nwave)
in.timebc = dt*3600*(cumsum(dum)-1);
%in.timebc = 100*(cumsum(dum));
%in.Tp= 4*dum;        % spectral peak period in seconds
in.Tp = T';
%in.Hrms = .3*dum;   % rms waveheightwave;  
in.Hrms = Hrms';

in.Wsetup = 0*dum;  % wave setup at seaward boundary in meters
in.swlbc = 0*dum;   % wave setup at seaward boundary in meters
in.angle = 0*dum;   % incident wave angle at seaward boundary in
in.jswl = 100; % Num of points below SWL (controls dx)

x0 = waves_short(1).x_wg(boundary_gauge);
dx = 1;
zb =[profile(1).z]-offset
x = [profile(1).x;]

x_1d =x0:dx:max(x);
zb_1d = interp1([0; x],[profile(1).z(1)-offset;zb],x_1d);
x = x-x0;
x_1d = x_1d-x0;
in.x = x_1d; % x points
in.zb = zb_1d % zb points
in.fw = .01*ones(size(in.zb)); % friction coeff
in.rwh = .01; 
in.wf = .03;
in.sg = 2.65;
in.d50 = .22;
in.effb = .005 % typical = .005 
in.efff = .01 % typical = .01 
in.slp = .2
in.tanphi=.66; %typical = .65
in.blp=.001; %typical .002


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%run cshore%%%%%%%%%%%%%%%%%%%%
if ~exist('cs_results')
  makeinfile(in);
  !/home/johnson/wes/cshore/code/cshore.out
  cs_results = load_results(in.nwave);
end
save results cs_results
