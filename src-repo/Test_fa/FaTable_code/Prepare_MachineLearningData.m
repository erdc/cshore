% clc
clear  

%%% wave and veg conditions
waterdepth  = 0.5 ;
H2h             = [0.005, 0.01:0.01:0.09, 0.1:0.05:0.75];
T                 = [0.6, 0.8, 1.0:0.5:4.0, 5.0:12.0] ; 
hv2h_total   = [0.1:0.1:0.8, 0.85:0.05:1.0, 1.1:0.1:1.5, 2.0] ; 
DATADIR    = './NewUr/' ; 

%%%               
Ursum        = [0.1:0.05:1, 1.5:0.5:1000]  ; 
mfinal       = zeros(size(Ursum))  ;
counthv2h    = 0 ; 
icolorcount  = 0 ; 
fid          = fopen('fa_database_headline.txt', 'w');
fprintf (fid, '##################################################\n') ;
fprintf (fid, 'hv/h       H/h        kh         Ur         fa\n') ; 

for ihv2h    = 1 : length(hv2h_total) 
    counthv2h    = counthv2h + 1 ;     
    hv2h         = hv2h_total (ihv2h) ; 
    
    modified_exp = zeros(size(Ursum))  ;
    count        = zeros(size(Ursum))  ;
    
    for iH2h = 1 : length (H2h)        
        waveheight = waterdepth *  H2h(iH2h) ;
        filename = [DATADIR, 'Oct4_h', strrep(num2str(waterdepth), '.', 'p'),...
                    '_H', strrep(num2str(waveheight), '.', 'p'),...
                    '_hv2h', strrep(num2str(hv2h), '.', 'p'),...
                    '_modified.txt'] ;
                
        data    = load (filename) ; 
        kh      = data (1, :) ;
        Ur      = data (2, :) ; 
        hv2h_exp_m  = data (3, :) ; 

        for iur = 1 : length (Ur)
            fprintf (fid, '%5.3f %10.6f %10.6f %10.6f %10.6f \n', ...
                          hv2h,   H2h(iH2h), kh(iur), Ur(iur), hv2h_exp_m(iur)) ; 
        end
    end
    
end

fclose (fid) ;
