clear 

fadatabase = load ('fa_database_noheadline.txt') ; 

%%%
hv2h_total   = [0.1:0.1:0.8, 0.85:0.05:1.0, 1.1:0.1:1.5, 2.0] ; 
H2h             = [0.005, 0.01:0.01:0.09, 0.1:0.05:0.75];
T                 = [0.6, 0.8, 1.0:0.5:4.0, 5.0:12.0] ; 
%%%

%%%
lhv2h = length (hv2h_total) ; 
lH2h  = length (H2h) ; 
lT       = length(T) ; 
%%%

testhv2h = 1.0246489287090943  ; 
testH2h = 0.19171997576854444   ;
testUr    = 11.456973451132177 ; 

%% find ihv2h
if testhv2h > hv2h_total(end)
    testhv2h = hv2h_total(end);
end

ihv2htmp = find ( hv2h_total >= testhv2h) ; % DONT USE MIN(ABS(array-element))
ihv2h_1 = ihv2htmp (1) 
if ihv2h_1 == 1
    ihv2h_1 = 2 ; 
end
ihv2h_0 = ihv2h_1 - 1 ;

fprintf ('hv/h:\n')
[hv2h_total(ihv2h_0), testhv2h, hv2h_total(ihv2h_1)]

%% find iH2h
iH2htmp = find ( H2h >= testH2h) ; % DONT USE MIN(ABS(array-element))iH2h_1 = iH2htmp (1) ;
iH2h_1 = iH2htmp (1) ;
if iH2h_1 == 1
    iH2h_1 = 2 ; 
end
iH2h_0 = iH2h_1 - 1  ;
% fprintf ('H/h:\n')
% [H2h(iH2h_0), testH2h, H2h(iH2h_1)]

%% Take out chunk of ihv2h_0 - ihv2h_1 and iH2h_0 - iH2h_1
% ihv2h_0, iH2h_0
chunk = fadatabase ((ihv2h_0-1)*lH2h*lT+1+(iH2h_0-1)*lT : (ihv2h_0-1)*lH2h*lT+1+(iH2h_0)*lT-1, :)  ;
chunkur = chunk (:, 4) ;
ichunkurtmp = find (chunkur>=testUr) ; 
ichunkur_1 = ichunkurtmp (1) ; 
ichunkur_0 = ichunkur_1 - 1 ; 
[chunkur(ichunkur_0), testUr, chunkur(ichunkur_1)]
% ihv2h_0, iH2h_0, ichunk_0 - ichunk_1
chunk_0001 = fadatabase ((ihv2h_0-1)*lH2h*lT+1+(iH2h_0-1)*lT+(ichunkur_0-1) : ...
                                            (ihv2h_0-1)*lH2h*lT+1+(iH2h_0-1)*lT+(ichunkur_0), [1,2,4,5])  ;
                                        
% ihv2h_0, iH2h_1
chunk = fadatabase ((ihv2h_0-1)*lH2h*lT+1+(iH2h_0)*lT : (ihv2h_0-1)*lH2h*lT+1+(iH2h_1)*lT-1, :) ;
chunkur = chunk (:, 4) ;
ichunkurtmp = find (chunkur>=testUr) ; 
ichunkur_1 = ichunkurtmp (1) ; 
ichunkur_0 = ichunkur_1 - 1 ; 
[chunkur(ichunkur_0), testUr, chunkur(ichunkur_1)]
% ihv2h_0, iH2h_1, ichunk_0 - ichunk_1
chunk_0101 = fadatabase ((ihv2h_0-1)*lH2h*lT+1+(iH2h_0)*lT+(ichunkur_0-1) : ...
                                            (ihv2h_0-1)*lH2h*lT+1+(iH2h_0)*lT+(ichunkur_0), [1,2,4,5])  ;

% ihv2h_1, iH2h_0
chunk = fadatabase ((ihv2h_0)*lH2h*lT+1+(iH2h_0-1)*lT : (ihv2h_0)*lH2h*lT+1+(iH2h_0)*lT-1, :)  ;
chunkur = chunk (:, 4) ;
ichunkurtmp = find (chunkur>=testUr) ; 
ichunkur_1 = ichunkurtmp (1) ; 
ichunkur_0 = ichunkur_1 - 1 ; 
% [chunk1ur(ichunkur_0), testUr, chunk1ur(ichunkur_1)]
% ihv2h_1, iH2h_0, ichunk_0 - ichunk_1
chunk_1001 = fadatabase ((ihv2h_0)*lH2h*lT+1+(iH2h_0-1)*lT+(ichunkur_0-1) : ...
                                            (ihv2h_0)*lH2h*lT+1+(iH2h_0-1)*lT+(ichunkur_0), [1,2,4,5])  ;
                                        
% ihv2h_1, iH2h_1
chunk = fadatabase ((ihv2h_0)*lH2h*lT+1+(iH2h_0)*lT : (ihv2h_0)*lH2h*lT+1+(iH2h_1)*lT-1, :)  ;
chunkur = chunk (:, 4) ;
ichunkurtmp = find (chunkur>=testUr) ; 
ichunkur_1 = ichunkurtmp (1) ; 
ichunkur_0 = ichunkur_1 - 1 ; 
% [chunk1ur(ichunkur_0), testUr, chunk1ur(ichunkur_1)]
% ihv2h_1, iH2h_0, ichunk_0 - ichunk_1
chunk_1101 = fadatabase ((ihv2h_0)*lH2h*lT+1+(iH2h_0)*lT+(ichunkur_0-1) : ...
                                            (ihv2h_0)*lH2h*lT+1+(iH2h_0)*lT+(ichunkur_0), [1,2,4,5])  ;

%% gather together and interpolate step by step
% [chunk_0001; chunk_0101; chunk_1001; chunk_1101]
fa00= interp1 (chunk_0001(:, 3), chunk_0001(:, 4), testUr, 'linear', 'extrap')  ;
fa01 = interp1 (chunk_0101(:, 3), chunk_0101(:, 4), testUr, 'linear', 'extrap')  ;
fa10 = interp1 (chunk_1001(:, 3), chunk_1001(:, 4), testUr, 'linear', 'extrap')  ;
fa11 = interp1 (chunk_1101(:, 3), chunk_1101(:, 4), testUr, 'linear', 'extrap')  ;

% chunk00_01 = [chunk_0001(1, 1:2), testUr, fa00;  ...
%                         chunk_0101(1, 1:2), testUr, fa01] 
% chunk10_11 = [chunk_1001(1, 1:2), testUr, fa10;  ...
%                         chunk_1101(1, 1:2), testUr, fa11] 
% 
% [hv2h_total(ihv2h_0), testhv2h, hv2h_total(ihv2h_1)]
% [H2h(iH2h_0), testH2h, H2h(iH2h_1)]

x0 = hv2h_total(ihv2h_0) ;
y0 = H2h(iH2h_0) ; 
x1 = hv2h_total(ihv2h_1) ;
y1 = H2h(iH2h_1) ;
fa_bilinearinterp = @(x, y) fa00 + (fa10-fa00) * ((x-x0)/(x1-x0)) + ...
                                                      (fa01-fa00) * ((y-y0)/(y1-y0)) + ...  
                                                      (fa11-fa01-fa10+fa00) * ((x-x0)/(x1-x0)) * ((y-y0)/(y1-y0)) ;
fa_test1 = fa_bilinearinterp (testhv2h, testH2h) ; 

[alpha] = SearchLookupTable (testH2h, testhv2h, testUr) ; 

[fa_test1, alpha]