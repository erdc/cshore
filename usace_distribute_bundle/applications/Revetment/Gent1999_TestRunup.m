% clc
clear all

%% This code adds veg effects on wave runup. this is paired with Gent_TestRunup.m


% global figID
addpath ../../../cshore/usace_distribute_bundle/mfiles
addpath ../../distribution/
addpath /Users/ling/Desktop/OngoingProjects/LingMatlabFuncs/
addpath ../Fv_SF_REQUIREDDATA

CaseType = 'B' ;

if strcmp(CaseType, 'A')
    Gent1999_ACases ;
    names = GentACases.name; 
    data = GentACases.data; 
elseif strcmp(CaseType, 'B')
    Gent1999_BCases ;
    names = GentBCases.name; 
    data = GentBCases.data; 
end

foreshoreslope = 1/100 ; 
dikeslope = 1/4 ; 
      
for icase = 1%:length(names) %7:9 %10:12 %%
    figID = icase ; 

    % script params
    iplotbc     = 0;                  % 1 to plot the applied boundary conditions
    iplotmorpho = 1;                  % 1 to plot the computed morphology results
    iplothydro  = 1;                  % 1 to plot the computed hydrodynamic results
    iplotrunup  = 1;                  % 1 to plot the computed runup position
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
    in.iweibull = 1;
    in.iveg   = 0;          % vegitation effect
                            % 0: no vegetation or vegetation represented by increased
                            % 1: veg. density, width, height and root depth are 
                            %    specified as input. The height and root depth 
                            %    vary with the bottom elevation change
                            % 2: veg. with constant density, width and height
    in.idiss  = 1;          % energy dissipation due to vegitation 
                            % (0: veg is accounted in bottom friction (original)
                            %  1: mendez, 2: chen_zhao)  
    in.veg_Cd = 2.0;        % vegitation drag coeff
    in.veg_Cdm = 2.0;        % vegitation drag coeff
    in.veg_n  = 420;       % vegitation density
    in.veg_dia= 0.008;    % vegitation diam
    in.veg_ht = 0.2;         % vegitation height
    in.veg_rod= 0.3;         % vegitation erosion limit below sand for failure
    in.veg_extent = [0.5425 0.86]; % vegitation coverage as fraction of total domain length
                            % [xveg_start, xveg_end] ./ total domain length
    in.dx     = 0.0075;       % constant dx 
    in.gamma  = 0.8;         % shallow water ratio of wave height to water depth
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
%     in.rwh    = .02;        % numerical rununp wire height 
    in.rwh    = 2.5e-3 ; 
    in.ilab   = 1;          % controls the boundary condition timing. Don't change

    % NOTE: if use in.ilab=0, error message shows up:
    % Data input is stopped because the start time for 
    % offshore wave conditions and water level is NOT ZERO

    in.fric_fac = .001;     % bottom friction factor
    
    
    % boundary conditions and timing
    ftime = 100;         % [sec] final time, dictates model duration
    dt = 120;            % time interval in seconds for wave and water level conditions
    if in.ilab==0 ;
      in.timebc_wave = [dt/2:dt:ftime];
    else
      in.timebc_wave = [0:dt:ftime];
    end
    in.timebc_surg = in.timebc_wave;
    in.nwave = length(in.timebc_wave); 
    in.nsurg = in.nwave; dum = ones(1,in.nwave);
    in.Tp= data(icase, 4)*dum ; % constant spectral peak period in seconds
%     in.Tp= Tps(icase)*dum ; % constant spectral peak period in seconds
% % %     in.Hrms = Hs1s(icase)/sqrt(2)*dum; % Hs=sqrt(2)*Hrms, Hrms=Hs/sqrt(2)
    in.Hrms = data(icase,3)/sqrt(2)*dum; % Hs=sqrt(2)*Hrms, Hrms=Hs/sqrt(2)
    in.Wsetup = 0*dum ;  % wave setup at seaward boundary in meters \bar{eta}
    in.swlbc = 0*dum;    % water level at seaward boundary in meters zb       
    
    in.angle = 0*dum;    % constant incident wave angle at seaward boundary in
    anguphase = 2*pi./in.Tp ; 
    in.freqmin= 0.1*anguphase ; 
    in.freqmax= 8.0*anguphase ; 
    in.numfreq= 500*ones(size(in.freqmax)) ; 

    % Idealized numerical tank
    Lx = 35.0;                  % length of domain
    flat_length = 30.0;         % length of flat portion at seaward end of numerical tank
    in.diketoe = 30.0 ; 
    zb_off = -(data(icase,1)-0.047);  % offshore bottom position (should be negative)
    zb_toe = zb_off + flat_length*foreshoreslope ; 
    zb_on = zb_toe + (Lx-flat_length)*dikeslope  ;            % onshore bottom position (should be pos)
    x = [0 flat_length Lx];     % x points
    zb = [zb_off zb_toe zb_on]; % zb points
    in.x = 0:in.dx:Lx;
    [j1 j2] = unique(x); 

    % original
    in.zb = interp1(x(j2),zb(j2),in.x);
%     figure
% plot (in.x, in.zb, '-k')
% return

    in.fw = in.fric_fac*ones(size(in.zb)); % cross-shore values of bot fric
    [kh] = dispersion(abs(zb_off), in.Tp(1))  ;
    wvlength = 2*pi/(kh/abs(zb_off)) ; 
    hv2h = in.veg_ht /abs(zb_off);

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%prepare input file for cshore%%%%%%%%%%%%%%%%%%%%
      makeinfile_usace_tryCdm(in) ;
      unix(['./../ICCE2018_Setup/sfapp']) ; 

%       makeinfile_usace_try(in) ;
%       unix(['./../ICCE2018/lwt']) ; 

%     makeinfile_usace(in);
%     unix(['./../origin']) ; 
%     unix(['./testMWL']) ; 
%     unix(['./nofvm']) ; 
    
    linecolor = 'b' ; 
    
    results = load_results_usace;

    if results.run_success==1;
      disp(' CSHORE computation completed ');
    elseif results.run_success==.5;
      disp(' CSHORE computation partially completed');
    end

    unix('rm -f O*') ; 
    if in.iveg==0
        z2pers_model_noveg(icase) = results.hydro.runup_2_percent ; 
        id = find (results.hydro.x<=30.0) ; 
% %     Hstoes_model_noveg(icase) = results.hydro.Hrms(id(end))*sqrt(2) ; 
        Hm0toes_model_noveg(icase) = results.hydro.sigma(id(end)) * 4 ;
    elseif in.iveg==1
        z2pers_model_veg(icase) = results.hydro.runup_2_percent ; 
    end

    if in.iveg==0
        figure(6); hold on; box on
%         subplot (2,1,1); hold on
%         plot (Hm0toes(icase), Hm0toes_model_noveg(icase), 'ok')
% 
%         subplot(2,1,2); hold on; box on
        plot (data(icase, 6), z2pers_model_noveg(icase), 'ob')
    end
end

if in.iveg==1
    fid=fopen(['Runup', CaseType, 'cases_iweibull', num2str(in.iweibull),  ...
               '_iveg', num2str(in.iveg), '.txt'], 'w') ;
    fprintf (fid, '%10.10f ', z2pers_model_veg) ;
    fclose (fid)
else
    fid=fopen(['Runup', CaseType, 'cases_iweibull', num2str(in.iweibull), ...
               '_iveg', num2str(in.iveg), '.txt'], 'w') ;
    fprintf (fid, '%10.10f ', z2pers_model_noveg) ;
    fclose (fid)
end

if in.iveg==0
    figure(6);hold on; box on
%     subplot(2,1,1); hold on; box on
%     xlabel ('measured H_{mo} (m)')
%     ylabel ('modeled H_{mo} (m)')
%     fplot(@(x) x, [0, 0.5], 'k')
%     fplot(@(x) 0.8*x, [0, 0.5], '--k')
%     fplot(@(x) 1.2*x, [0, 0.5], '--k')
%     xlim([0, 0.2])
%     ylim([0, 0.2])
%     set(gca, 'fontsize', 12)
% 
%     subplot(2,1,2); hold on; box on
    xlabel ('measured R_{2%} (m)')
    ylabel ('modeled R_{2%} (m)')
    fplot(@(x) x, [0, 0.6], 'k')
    fplot(@(x) 0.8*x, [0, 0.6], '--k')
    fplot(@(x) 1.2*x, [0, 0.6], '--k')
    xlim([0, 0.6])
    ylim([0, 0.6])
    set(gca, 'fontsize', 12)    
end

return
figure; 
subplot(2,1,1); hold on; box on
% plot (z2pers_model_noveg, z2pers_model_veg, 'or')
xlabel ('modeled R_{2%} without veg (m)')
ylabel ('modeled R_{2%} with veg (m)')
fplot(@(x) x, [0, 0.6], 'k')
fplot(@(x) 0.8*x, [0, 0.6], '--k')
fplot(@(x) 0.5*x, [0, 0.6], '--k')
xlim([0, 0.6])
ylim([0, 0.6])
set(gca, 'fontsize', 12)


rmpath ../../cshore/usace_distribute_bundle/mfiles
rmpath ../../distribution/
rmpath /Users/ling/Desktop/OngoingProjects/LingMatlabFuncs/

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
