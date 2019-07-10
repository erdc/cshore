% clc
clear 

% global figID
addpath ../../mfiles
addpath ./USDA_data/wave_setup

% load lab data
DIR = ['./USDA_data/wave_setup/']; 

cases=dir([DIR, 'RANDOM*.dat']) ; 

% offshore w/ veg
mwl_offshore_veg = [0.25*1e-3, 0.1*1e-3, 0*1e-3, 0*1e-3, -.5*1e-3, -.5*1e-3, -.3*1e-3] ; 

% offshore w/o veg
% mwl_offshore_noveg = [0.25*1e-3, 0.1*1e-3, 0*1e-3, 0*1e-3, -.5*1e-3, -.5*1e-3, -.3*1e-3] ; 
mwl_offshore_noveg = [0.25*1e-3, 0.1*1e-3, 0*1e-3, 0*1e-3, -.5*1e-3, -.5e-3, -.3*1e-3] ; 

CDs  =   [2.3, 2.0,  2.0, 2.0,  1.7,  1.6,  1.6]; 
CDMs = [1.9, 2.0,  1.6, 2.0,  1.2,  1.5,  0.9]; 
khs  = [1.3, 0.89, 1.3, 0.77, 0.89, 0.55, 0.77];

for icase = 1 : length(cases) 
    
    filename = cases(icase).name 
    
    text   = fileread([DIR, filename]) ; 
    id1 = strfind (text, '(s)') ; 
    id2 = strfind (text, 'Water Depth') ;
    PeakPer = text(id1+6:id2-1) ; 
    PeakPer = strrep (PeakPer, ' ', '') ;
    
    fid = fopen ([DIR, filename], 'rt')  ;
    labdata = textscan(fid,'%f%f%f%f%f%f%f%f%f', ...
            'Delimiter',' ','MultipleDelimsAsOne',true,'Headerlines',5);
    figID = icase ; 

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
    in.iweibull = 0 ;
    in.iveg   = 3;          % vegitation effect
                            % 0: no vegetation or vegetation represented by increased
                            % 1: veg. density, width, height and root depth are 
                            %    specified as input. The height and root depth 
                            %    vary with the bottom elevation change
                            % 2: veg. with constant density, width and height
    in.veg_Cd = CDs(icase);        % vegitation drag coeff   
    in.veg_Cdm = in.veg_Cd ;    % for iFv ==2, Cdm is specified by CDMs                            
if in.iveg==3                        
    in.idiss  = 1;          % energy dissipation due to vegitation 
                            % (0: veg is accounted in bottom friction (original)
                            %  1: mendez, 2: chen_zhao)  
    in.iFv = 2 ;
    if (in.iFv == 2)
        in.veg_Cdm = CDMs(icase);    
    end
    
end
    
    in.veg_n  = 3150;       % vegitation density
    in.veg_dia= 3.175e-3;       % vegitation diam
    in.veg_ht = 0.2;         % vegitation height
    in.veg_rod= 0.3;         % vegitation erosion limit below sand for failure
    in.veg_extent = [0.5425 0.7151]; % vegitation coverage as fraction of total domain length
                            % [xveg_start, xveg_end] ./ total domain length
    in.dx     = 0.005;       % constant dx 
    in.gamma  = .7;         % shallow water ratio of wave height to water depth
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

    in.fric_fac = .02;     % bottom friction factor
    
    
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
    in.Tp= str2num(PeakPer)*dum ; % constant spectral peak period in seconds
    if (in.iveg==0)
        in.Hrms = labdata{4}(1)/sqrt(2)*dum; % Hs=sqrt(2)*Hrms, Hrms=Hs/sqrt(2)
%         in.Wsetup = labdata{2}(1)*1e-3*dum;   % wave setup at seaward boundary in meters
%         in.swlbc = labdata{2}(1)*1e-3*dum;    % water level at seaward boundary in meters
        in.Wsetup = mwl_offshore_noveg(icase)*dum ; %%labdata{6}(2)*1e-3*dum ;  % wave setup at seaward boundary in meters \bar{eta}
        in.swlbc = 0*dum;    % water level at seaward boundary in meters zb       
    else
        in.Hrms = labdata{8}(1)/sqrt(2)*dum; % Hs=sqrt(2)*Hrms, Hrms=Hs/sqrt(2)
        in.Wsetup = mwl_offshore_veg(icase)*dum ; %%labdata{6}(2)*1e-3*dum ;  % wave setup at seaward boundary in meters \bar{eta}
        in.swlbc = 0*dum;    % water level at seaward boundary in meters zb       
    end
    
    in.angle = 0*dum;    % constant incident wave angle at seaward boundary in
    anguphase = 2*pi./in.Tp ; 
    if in.iveg==3   && in.idiss==2 
        in.freqmin= 0.1*anguphase ; 
        in.freqmax= 8.0*anguphase ; 
        in.numfreq= 500*ones(size(in.freqmax)) ; 
    end
    
    % Idealized numerical tank
    Lx = 21.2;                  % length of domain
    zb_off = -0.4;              % offshore bottom position (should be negative)
    zb_on = 0.2667;                 % onshore bottom position (should be pos)
    flat_length = 7.2;         % length of flat portion at seaward end of numerical tank
    x = [0 flat_length Lx];     % x points
    zb = [zb_off zb_off zb_on]; % zb points
    in.x = 0:in.dx:Lx;
    [j1 j2] = unique(x); 

    % original
    in.zb = interp1(x(j2),zb(j2),in.x);
    
    % modify  % inscribed circle of a triangle
    id    = find(in.x>flat_length); 
    x0    = in.x(id(1))   ;  
    nnx1  = 200 ;
    nnx2  = 200 ; 
    x1    = in.x(id(1)-nnx1) ; 
    x2    = in.x(id(1)+nnx2) ; 
    y2    = (x2-x0)/21.0 ; 
    y1    = ((x2-x1)^2 + y2^2)/2/y2 ; 
    circ  = @(x) zb_off+y1 - sqrt(y1^2 - (x-x1).^2) ;
    
    id2   = id(1)-nnx1+1:id(1)+nnx2-1 ; 
    zb1   = in.zb (1:id(1)-nnx1) ;
    zb3   = in.zb (id(1)+nnx2:length(in.x)) ;
    
    zb2   = circ(in.x(id2)) ; 
    in.zb = [zb1, zb2, zb3] ; 
      
    in.fw = in.fric_fac*ones(size(in.zb)); % cross-shore values of bot fric

    [kh] = dispersion(abs(zb_off), in.Tp(1))  ;
    wvlength = 2*pi/(kh/abs(zb_off)) ; 
    hv2h = in.veg_ht /abs(zb_off);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%prepare input file for cshore%%%%%%%%%%%%%%%%%%%%
    makeinfile_usace_vegfeature(in) ;
    unix(['./../../../src-repo/updatedveg']) ;  
    
    results = load_results_usace;
    if results.run_success==1;
      disp(' CSHORE computation completed ');
    elseif results.run_success==.5;
      disp(' CSHORE computation partially completed');
    end

%     plot_results_usace
    figure(figID+200);

    subplot (2,1,1); hold on; box on
    set (gca, 'fontsize', 15)
    if (in.iveg==0)
        errorbar (labdata{1}, labdata{4}/sqrt(2)./in.Hrms, labdata{5}/sqrt(2)./in.Hrms, 'ok')    
    else
        errorbar (labdata{1}, labdata{8}/sqrt(2)./in.Hrms, labdata{9}/sqrt(2)./in.Hrms, 'ok')    
    end
    plot (results.hydro.x,results.hydro.Hrms./in.Hrms, 'k', 'linewidth',2) 
%     legend ('measurements', 'modeled')
    if (in.iveg==1)
        plot (in.veg_extent(1)*Lx*ones(1,2), linspace(0, 1.5, 2), ':k', 'linewidth', 2)
        plot (in.veg_extent(2)*Lx*ones(1,2), linspace(0, 1.5, 2), ':k', 'linewidth', 2)
    end
    xlim ([2.5, 15.3])
    ylim ([0, 1.5])
    ylabel ('H_{rms}./H_{rms, 0} ')
    xlabel ('x (m)')
    id=find(in.x>in.veg_extent(1)); 
    H2hs(icase) = results.hydro.Hrms(id(1))/results.hydro.depth(id(1)) ;
    
    casename = strrep(strrep(strrep(filename, '.dat', ''), '_', ' '), 'RANDOM ', '') ;
    title (['Case: ', casename])
    
    subplot (2,1,2); hold on; box on
    set (gca, 'fontsize', 15)
    if (in.iveg==0)
        errorbar (labdata{1}, labdata{2}*0.001./in.Hrms, labdata{3}*0.001./in.Hrms, 'ok')    
    else
        errorbar (labdata{1}, labdata{6}*0.001./in.Hrms, labdata{7}*0.001./in.Hrms, 'ok')    
    end
    plot (results.hydro.x, results.hydro.setup./in.Hrms, '-k', 'linewidth', 2)
%     legend ('measurements', 'modeled')
    if (in.iveg==1)
        plot (in.veg_extent(1)*Lx*ones(1,2), linspace(-0.1, 0.5, 2), ':k', 'linewidth', 2)
        plot (in.veg_extent(2)*Lx*ones(1,2), linspace(-0.1, 0.5, 2), ':k', 'linewidth', 2)
    end
    xlim ([2.5, 15.3])
    ylabel ('$$\overline{\eta} /H_{rms,0}$$','interpreter','latex')
    xlabel ('x (m)')
    ylim ([-0.1, 0.05])
    
    %% calibrate setup
    idx     = find (labdata{1}<in.veg_extent(1)*Lx) ; 
    idx     = idx (2:end-1) ; 
    exp_x   = labdata{1}(idx) ; 
    exp_mwl = labdata{6}(idx) ; 
    
    for ii = 1:length(idx)
        id1     = find(results.hydro.x<exp_x(ii)) ;
        num_x(ii) = results.hydro.x(id1(end)) ;
        num_mwl(ii) = results.hydro.setup(id1(end))*1000 ; 
    end
        
    error = sqrt(sum((num_mwl - exp_mwl').^2) ./ sum(exp_mwl.^2))  
    
    %% Plot MWL_measure VS MWL_model
    idx1     = find (labdata{1}>=in.veg_extent(1)*Lx & labdata{1}<=in.veg_extent(2)*Lx) ; 
    xgage    = labdata{1}(idx1) ;
    MWL_lab  = labdata{6}(idx1) ; 
    MWL_lab_STD = labdata{7}(idx1) ; 
    
    MWL_model = zeros(size(MWL_lab)) ; 
    for ii = 1:length(idx1)
        [val1, id1]  = min(abs(results.hydro.x-xgage(ii))) ;
        MWL_model(ii)    = results.hydro.setup(id1) * 1000 ; 
    end
        
    figure(1124); hold on; box on
    hh = herrorbar (MWL_lab, MWL_model, -MWL_lab_STD, MWL_lab_STD, 'ok')  ; 
    axis square
    
%     fprintf (fid, '%10.10e ', MWL_model, MWL_model)
    
    %% Plot H_measure VS H_model
    idx1     = find (labdata{1}>=in.veg_extent(1)*Lx & labdata{1}<=in.veg_extent(2)*Lx) ; 
    xgage    = labdata{1}(idx1) ;
    H_lab    = labdata{8}(idx1)/sqrt(2)*100 ;
    H_lab_STD= labdata{9}(idx1)/sqrt(2)*100 ;
    
    H_model = zeros(size(H_lab)) ; 
    for ii = 1:length(idx1)
        [val1, id1]  = min(abs(results.hydro.x-xgage(ii))) ;
        H_model(ii)  = results.hydro.Hrms(id1)*100 ; 
    end
    
    figure(1125); hold on; box on
    hh = herrorbar (H_model, H_lab, -H_lab_STD, H_lab_STD, 'ok')  ; 
    axis square

    unix('rm -f O*') ; 
%     unix('rm -f *.DAT') ; 

end

figure(1124);
fplot (@(x) x, 'k', 'linewidth', 1)
xlabel ('measured MWL (mm)')
ylabel ('modeled MWL (mm)')
xlim ([-3, 4])
ylim ([-3, 4])
set (gca, 'fontsize', 15)
axis square

figure(1125);
fplot (@(x) x, 'k', 'linewidth', 1)
xlabel ('measured H_{rms} (cm)')
ylabel ('modeled H_{rms} (cm)')
xlim ([0, 8])
ylim ([0, 8])
set (gca, 'fontsize', 15)
axis square


rmpath ../../mfiles
rmpath ./USDA_data/wave_setup

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
