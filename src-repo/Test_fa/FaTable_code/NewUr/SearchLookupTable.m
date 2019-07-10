function [alpha] = SearchLookupTable (Hs2h, hv2h, Ur)
% [fa] = SearchLookupTable (H2h, hv2h, Ur)

    hv2h_table  = [0.1:0.1:0.8, 0.85:0.05:1.0, 1.1:0.1:1.5, 2.0] ; 
    H2h_table   = [0.005, 0.01:0.01:0.09, 0.1:0.05:0.75];
    T_table        = [0.6, 0.8, 1.0:0.5:4.0, 5.0:12.0] ; 
    DATADIR     = ['/Users/lzhu/Desktop/OngoingProjects/4_WaveSetupDragVeg/Codes/TUDenmark_code/NewUr/'] ; 

    %%% search lookup table for alpha
    if Hs2h <= H2h_table(1)
        id_H2h = 1  ;
    else
        tmp         = find (H2h_table<=Hs2h) ; 
        id_H2h      = tmp (end) ; 
    end
    
    tmp         = find (hv2h_table<=hv2h)  ;
    id_hv2h     = tmp (end) ; 
    if hv2h<hv2h_table(1)
        id_hv2h = 1 ; 
    end
    
    if id_H2h == 1 || id_H2h==length(H2h_table)
        waveheight = H2h_table (id_H2h) * 0.5 ;
        if id_hv2h == 1 || id_hv2h == length(hv2h_table)
            'case1' ;
            filename1  = [DATADIR, 'Oct4_h', strrep(num2str(0.5), '.', 'p'),...
                '_H', strrep(num2str(waveheight), '.', 'p'),...
                '_hv2h', strrep(num2str(hv2h_table (id_hv2h)), '.', 'p'),...
                '_modified.txt']  ;
            data1      = load (filename1) ; 
            Ur1        = data1(2, :) ;
            alpha1     = data1(3, :) ;
            [Ur1, index] = unique(Ur1); 
            alpha      = interp1(Ur1, alpha1(index), Ur) ;                
        else
            'case2'  ;
            filename1  = [DATADIR, 'Oct4_h', strrep(num2str(0.5), '.', 'p'),...
                '_H', strrep(num2str(waveheight), '.', 'p'),...
                '_hv2h', strrep(num2str(hv2h_table (id_hv2h)), '.', 'p'),...
                '_modified.txt']  ;
            data1      = load (filename1) ; 
            Ur1        = data1(2, :) ;
            alpha1     = data1(3, :) ;
            [Ur1, index] = unique(Ur1); 
            alpha_val1 = interp1(Ur1, alpha1(index), Ur) ; 
            filename2  = [DATADIR, 'Oct4_h', strrep(num2str(0.5), '.', 'p'),...
                '_H', strrep(num2str(waveheight), '.', 'p'),...
                '_hv2h', strrep(num2str(hv2h_table (id_hv2h+1)), '.', 'p'),...
                '_modified.txt']  ;
            data2      = load (filename2) ; 
            Ur2        = data1(2, :) ;
            alpha2     = data1(3, :) ;
            [Ur2, index] = unique(Ur2); 
            alpha_val2 = interp1(Ur2, alpha2(index), Ur) ;                
            alpha      = interp1(hv2h_table(id_hv2h:id_hv2h+1), [alpha_val1, alpha_val2], hv2h) ;            
        end 
    else
        if id_hv2h == 1 || id_hv2h == length(hv2h_table)
            'case3' ;
            waveheight = H2h_table (id_H2h) * 0.5 ;
            filename1  = [DATADIR, 'Oct4_h', strrep(num2str(0.5), '.', 'p'),...
                '_H', strrep(num2str(waveheight), '.', 'p'),...
                '_hv2h', strrep(num2str(hv2h_table (id_hv2h)), '.', 'p'),...
                '_modified.txt']  ;
            waveheight = H2h_table (id_H2h+1) * 0.5 ;
            filename2  = [DATADIR, 'Oct4_h', strrep(num2str(0.5), '.', 'p'),...
                '_H', strrep(num2str(waveheight), '.', 'p'),...
                '_hv2h', strrep(num2str(hv2h_table (id_hv2h)), '.', 'p'),...
                '_modified.txt']  ;
            data1      = load (filename1) ; 
            Ur1        = data1(2, :) ;
            alpha1     = data1(3, :) ;   
            [Ur1, index] = unique(Ur1); 
            alpha_val1 = interp1(Ur1, alpha1(index), Ur) ;  

            data2      = load (filename2) ; 
            Ur2        = data2(2, :) ;
            alpha2     = data2(3, :) ;
            [Ur2, index] = unique(Ur2); 
            alpha_val2 = interp1(Ur2, alpha2(index), Ur) ;
            alpha      = interp1(H2h_table (id_H2h:id_H2h+1), [alpha_val1, alpha_val2], Hs2h, 'linear', 'extrap') ;            
            
        else
            'case4' ;
            waveheight = H2h_table (id_H2h) * 0.5 ;
            filename1  = [DATADIR, 'Oct4_h', strrep(num2str(0.5), '.', 'p'),...
                '_H', strrep(num2str(waveheight), '.', 'p'),...
                '_hv2h', strrep(num2str(hv2h_table (id_hv2h)), '.', 'p'),...
                '_modified.txt']  ;
            data1      = load (filename1) ; 
            Ur1        = data1(2, :) ;
            alpha1     = data1(3, :) ;
            [Ur1, index] = unique(Ur1);                                 
            alpha_val1 = interp1(Ur1, alpha1(index), Ur) ;

            filename2  = [DATADIR, 'Oct4_h', strrep(num2str(0.5), '.', 'p'),...
                '_H', strrep(num2str(waveheight), '.', 'p'),...
                '_hv2h', strrep(num2str(hv2h_table (id_hv2h+1)), '.', 'p'),...
                '_modified.txt'] ;
            data2      = load (filename2) ; 
            Ur2        = data2(2, :) ;
            alpha2     = data2(3, :) ;
            [Ur2, index] = unique(Ur2);                 
            alpha_val2 = interp1(Ur2, alpha2(index), Ur) ;

            waveheight = H2h_table (id_H2h+1) * 0.5 ;  
            filename3  = [DATADIR, 'Oct4_h', strrep(num2str(0.5), '.', 'p'),...
                '_H', strrep(num2str(waveheight), '.', 'p'),...
                '_hv2h', strrep(num2str(hv2h_table (id_hv2h)), '.', 'p'),...
                '_modified.txt'] ;
            data3      = load (filename3) ; 
            Ur3        = data3(2, :) ;
            alpha3     = data3(3, :) ;
            [Ur3, index] = unique(Ur3); 
            alpha_val3 = interp1(Ur3, alpha3(index), Ur) ;

            filename4  = [DATADIR, 'Oct4_h', strrep(num2str(0.5), '.', 'p'),...
                '_H', strrep(num2str(waveheight), '.', 'p'),...
                '_hv2h', strrep(num2str(hv2h_table (id_hv2h+1)), '.', 'p'),...
                '_modified.txt'] ;
            data4      = load (filename4) ; 
            Ur4        = data4(2, :) ;
            alpha4     = data4(3, :) ;  
            [Ur4, index] = unique(Ur4);                 
            alpha_val4 = interp1(Ur4, alpha4(index), Ur) ;
            valtmp     = [alpha_val1, alpha_val3; alpha_val2 , alpha_val4] ; 
            alpha      = interp2(H2h_table(id_H2h:id_H2h+1), hv2h_table(id_hv2h:id_hv2h+1), valtmp, Hs2h, hv2h) ;   
        end
    end    

end