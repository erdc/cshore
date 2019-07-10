clear 

fdir='/Users/lzhu/Desktop/OngoingProjects/1_EMI_MomentumBalance/Results/PureWave_T1/scpout/';
outputdir = '/Users/lzhu/Desktop/OngoingProjects/1_EMI_MomentumBalance/Results/PureWave_T1/' ; 

m=2048;
dx=0.1 ; 
xgrid_all =[0:m-1]*dx; % bathy
NZ         = 20 ; 

% find out xgrid ids for study area
studyextend = [10.0, 190.0] ; % in meter 
studyid = find (xgrid_all>studyextend(1) & xgrid_all<studyextend(2)) ; 
% NOTE: studyid must be continuous, otherwise, cannot perform derivative.
xgrid = xgrid_all (studyid) ;

% find out xgrid for vegetation
vegextend = [130.0, 180.0] ; % in meter 
vegid = find (xgrid_all>vegextend(1) & xgrid_all<vegextend(2)) ; 

% put wave period
waveper = 3.2 ;  % in sec
omega   = 2*pi / waveper ;
rho         = 1000.0 ; 
grav       = 9.806 ;

% load depth 
depthtmp = load ([fdir, 'depth.txt']) ;
depth       = depthtmp (2, studyid) ;

% load time
time_all = load ([fdir, 'time_dt.out']) ; 
timeid = 2000:8001 ;
time       = time_all (timeid, 1) ; 
NT         = floor ( ( time (end) - time (1) ) / waveper ) ; 

%2nd order accurate difference matrices. 
D1 = fdmatrix(xgrid, 1, 2);

%
clearvars eta_all
eta_all = []; 
for itime = timeid
    
    etatmp   = load ([fdir, 'eta_', num2str(itime, '%05.f')]) ; 
    eta         = etatmp (2, studyid)  ; 
    eta_all   = [eta_all; eta] ; 
end

% Separate incident and reflected waves and store them to files.
[Sin,Sout]=RadonSeparation(M) ;
save ([outputdir, 'Sin.mat'], 'Sin')
save ([outputdir, 'eta_all.mat'], 'eta_all')

%% compute MWL 
% 1. get integer folds of wave period
t1           = time (1) ;
t2           = time (1) + NT*waveper ; 
[~, t2idtmp]  = min ( abs(time - t2) ) ;
t2id        = t2idtmp(end) ; 

% 2. compute MWL 
for ix = 1 : length (xgrid) 
    etabar (ix) = mean (eta_all(1:t2id, ix)) ; %trapz (time(1:t2id), Sxx (1:t2id, ix)) / (time(t2id) - time(1)) ;
    etabar_2 (ix) = trapz (time(1:t2id), eta_all (1:t2id, ix)) / (time(t2id) - time(1)) ;
end

%% compute MWL changes
detabardx = D1 * etabar' ;
LHS = rho*grav*(depth+etabar).*detabardx' ;

%% Retrieve time series of eta at x = 120 m and 130 m
id =find (xgrid >= 100 & xgrid <=150) ; 

%% Plot 
plot (xgrid, etabar', '-k', 'linewidth', 2)
plot (xgrid, etabar_2', '-c', 'linewidth', 2)
legend ('$\overline{\eta}$', '$\overline{\eta} \; 2$', 'Interpreter', 'latex')
xlim([0, 180])
set (gca, 'fontsize', 15)
xlabel ('x (m)')
ylabel ('$\overline{\eta}\;\mbox{(m)}$', 'Interpreter', 'latex')

fid = fopen ([outputdir, 'MWL.txt'], 'w') ;
fprintf (fid, '%10.10f ', xgrid) ;
fprintf (fid, '\n') ; 
fprintf (fid, '%10.10f ', etabar) ;
fclose (fid) ;

fid = fopen ([outputdir, 'MWL_LHS.txt'], 'w') ;
fprintf (fid, '%10.10f ', xgrid) ;
fprintf (fid, '\n') ; 
fprintf (fid, '%10.10f ', LHS) ;
fclose (fid) ;

fid = fopen ([outputdir, 'eta_x100m.txt'], 'w') ;
fprintf (fid, '%10.10f ', time) ; 
fprintf (fid, '\n') ; 
fprintf (fid, '%10.10f ', eta_all(:, id(1))) ; 
fclose (fid) ;

fid = fopen ([outputdir, 'eta_x150m.txt'], 'w') ;
fprintf (fid, '%10.10f ', time) ; 
fprintf (fid, '\n') ; 
fprintf (fid, '%10.10f ', eta_all(:, id(end))) ; 
fclose (fid) ;


