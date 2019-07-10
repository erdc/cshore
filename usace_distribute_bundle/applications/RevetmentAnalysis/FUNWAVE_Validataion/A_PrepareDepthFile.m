clear 

vegtest = 1 ;

% load lab data
revetment1 = struct ('id', 'RS20B1', 'd1', 38.9/100, 'Tp',  2.3, 'Hrms', 11.2/100, 'dt', 20.6/100, 'Rc', 18.6/100, 'etabar', 0.28/100, 'SWL', 0, 'R2pmeasure', 13.31/100) ; 
revetment2 = struct ('id', 'RS20C1', 'd1', 38.9/100, 'Tp',  3.0, 'Hrms', 7.3/100, 'dt', 20.6/100, 'Rc', 18.6/100, 'etabar', 0.17/100, 'SWL', 0, 'R2pmeasure', 13.59/100) ; 
revetment3 = struct ('id', 'RS22B1', 'd1', 40.9/100, 'Tp',  2.3, 'Hrms', 11.6/100, 'dt', 22.6/100, 'Rc', 16.6/100, 'etabar', 0.11/100, 'SWL', 2/100,  'R2pmeasure', 14.07/100) ; 
revetment4 = struct ('id', 'RS22C1', 'd1', 40.9/100, 'Tp',  2.9, 'Hrms', 7.6/100, 'dt', 22.6/100, 'Rc', 16.6/100, 'etabar', 0.13/100, 'SWL', 2/100,  'R2pmeasure', 14.60/100) ; 
revetment5 = struct ('id', 'RS24B1', 'd1', 42.9/100, 'Tp',  2.3, 'Hrms', 11.9/100, 'dt', 24.6/100, 'Rc', 14.6/100, 'etabar', -0.27/100, 'SWL', 4/100,  'R2pmeasure', 14.12/100) ; 
revetment6 = struct ('id', 'RS24C1', 'd1', 42.9/100, 'Tp',  2.9, 'Hrms', 7.8/100, 'dt', 24.6/100, 'Rc', 14.6/100, 'etabar', 0, 'SWL', 4/100, 'R2pmeasure', 14.98/100) ;   

for icase = 6 %1 : 6
    eval (['in.Hrms = revetment', num2str(icase), '.Hrms;'] )   ; % Hs=sqrt(2)*Hrms, Hrms=Hs/sqrt(2)
    eval (['in.Tp = revetment', num2str(icase), '.Tp ;']) ;  % constant spectral peak period in seconds
    eval (['d1 = revetment', num2str(icase), '.d1 ;']) ; 
    eval (['dt = revetment', num2str(icase), '.dt ;']) ; 
    eval (['Rc = revetment', num2str(icase), '.Rc ;']) ; 
    
    %% Prepare Depth File
     x_p         = [0,                     6.3,      6.3+0.82,                                6.3+0.82+1.24] ;
%     zb_p     = [-dt-0.1831,        -dt,       -dt+0.02384,                            Rc-0.119] ; 
    zb_p     = [-dt-0.1831,        -dt,       -dt+0.02384,                            Rc-0.25] ; 
    
    % Extend horizontal parts before x_p
    x_p  = x_p + 40 ;
    x_p = [0, x_p] ; 
    zb_p = [zb_p(1), zb_p] ; 

    Lx       = 60;  % length of domain
    in.dx     = 0.005;       % constant dx 
    in.x_p   = 0:in.dx:Lx;
    in.x       = 0:in.dx:Lx;
    in.zb_p = interp1(x_p, zb_p, in.x_p, 'linear', 'extrap');  
    
    in.zb     = in.zb_p; 
        
    [kh] = dispersion (d1, in.Tp) ;
    k = kh / d1 ; L = 2*pi/k ;
    M = 2048 ; 
    N = 3 ;
    DX = 0.03 ;
    DY = 1.0 ;
    Xgrid = 0:DX:(M-1)*DX ; 
    Depth_tmp = interp1(in.x, in.zb, Xgrid, 'linear', 'extrap') ; 
    Depth = repmat (Depth_tmp, [N, 1]) ; 
    % For FUNWAVE-TVD, depth in water should be positive
    Depth = - Depth ;
    fid0 = fopen (['depth_case', num2str(icase), '.txt'], 'w') ; 
    for ii = 1 : N
        for jj = 1:M
            fprintf (fid0, '%5.5f ', Depth(ii, jj)) ;
        end
        fprintf (fid0, '\n') ;
    end
    fclose (fid0);

    fid00 = fopen (['xgrid_case', num2str(icase), '.txt'], 'w') ; 
    fprintf (fid00, '%5.5f ', Xgrid) ;
    fclose (fid00);
    
    
    %% Prepare Vegetation Files
    if vegtest ==1
        vegextend = [6.3, (6.3+0.82+1.24)] + 40 ;
    else 
        vegextend = [1.96, (6.3+0.82+1.24)] + 40 ;
    end 

    xid = find (Xgrid>vegextend(1) & Xgrid<vegextend(2)) ; 

    vegheight = zeros(M, N) ; 
    vegNv = zeros(M, N) ; 
    vegbv = zeros(M, N) ; 
    vegCd = zeros(M, N) ; 
    
    vegheight(xid, :) = 0.2 ; 
    vegNv(xid, :) = 400.0; 
    vegbv(xid, :) = 0.01; 
    vegCd(xid, :) = 1.8578; 

    fid1 = fopen (['veghv_case', num2str(icase), '.txt'], 'w') ; 
    fid2 = fopen (['vegNv_case', num2str(icase), '.txt'], 'w') ; 
    fid3 = fopen (['vegbv_case', num2str(icase), '.txt'], 'w') ; 
    fid4 = fopen (['vegCd_case', num2str(icase), '.txt'], 'w') ; 
    for ii = 1:N
        for jj = 1:M
            fprintf (fid1, '%5.5f ', vegheight(jj, ii)) ;
            fprintf (fid2, '%5.5f ', vegNv(jj, ii));
            fprintf (fid3, '%5.5f ', vegbv(jj, ii));
            fprintf (fid4, '%5.5f ', vegCd(jj, ii));
        end
        fprintf (fid1, '\n') ;
        fprintf (fid2, '\n') ;
        fprintf (fid3, '\n') ;
        fprintf (fid4, '\n') ;
    end
    fclose (fid1);
    fclose (fid2);
    fclose (fid3);
    fclose (fid4);    
end

