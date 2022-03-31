function makeinfile_usace_flexveg_CdEPH2CSHORE(in)
% function makeinfile_usace(in)

%First check for NaN
if max(isnan([in.x(:);in.zb(:);in.Tp(:);in.Hrms(:);in.swlbc(:);in.Wsetup(:);in.angle(:)]))
  error('The in structure contains NaN')
end

% if in.iveg~=3, in.ivegCd and and in.ivegtype have to be 0
if in.iveg~=3
    in.ivegCd = 0 ;
    in.ivegtype = 0 ;
end

fid = fopen('infile','w');
fprintf(fid,'%i \n',length(in.header));
for i = 1:length(in.header)
  fprintf(fid,'%s \n',cell2mat(in.header(i)));
end
fprintf(fid,'%-8i                                  ->ILINE\n',in.iline);
fprintf(fid,'%-8i                                  ->IPROFL\n',in.iprofl);
if floor(in.iprofl)==1
    fprintf(fid,'%-8i                                  ->ISEDAV\n',in.isedav);
end
% %note:
% C     IPROFL=0 for fixed bottom profile(no input for ISEDAV=0)
% C     IPROFL=1 for profile evolution computation(input ISEDAV)
% C     IPROFL=1.1 for IPROFL=1 and no initial bottom smoothing
% C     IPROFL=2 for dike erosion computation (ISEDAV=0)
if floor(in.iprofl)==0
    in.isedav = 0 ;
end

fprintf(fid,'%-8i                                  ->IPERM\n',in.iperm);
fprintf(fid,'%-8i                                  ->IOVER\n',in.iover);
if in.iover
  fprintf(fid,'%-8i                                  ->IWTRAN\n',in.iwtran);
  if in.iwtran==0
    fprintf(fid,'%-8i                                  ->IPOND\n',in.ipond);
  end
end
if in.iover==1&in.iperm==0&floor(in.iprofl)==1
  fprintf(fid,'%-8i                                  ->INFILT\n',in.infilt);
end
fprintf(fid,'%-8i                                  ->IWCINT\n',in.iwcint);
fprintf(fid,'%-8i                                  ->IROLL \n',in.iroll);
fprintf(fid,'%-8i                                  ->IWIND \n',in.iwind);
fprintf(fid,'%-8i                                  ->ITIDE \n',in.itide);
fprintf(fid,'%-8i                                  ->IVEG  \n',in.iveg);

if in.iveg==3
   fprintf(fid,'%-8i                                  ->IDVEGCD  \n',in.ivegCd);
   fprintf(fid,'%-8i                                  ->IDVEGTYPE  \n',in.ivegtype);

   if isfield(in, 'Cdcap')
       fprintf(fid,'%11.4f                              ->CdCap  \n',in.Cdcap);
   else %%
       fprintf(fid,'%11.4f                              ->CdCap  \n',45.0);
   end

   if isfield(in, 'rhowater')
       fprintf(fid,'%11.4f                              ->rhowater  \n',in.rhowater);
   else %% 
       fprintf(fid,'%11.4f                              ->rhowater  \n',1000.0);
   end
   
% note: use "IDVEGCD" and "IDVEGTYPE" instead of "IVEGCD" and "IVEGTYPE" so that in 
% load_results_usace, it will only read IVEG line.

   if isfield(in, 'nvegsegment')
       fprintf(fid,'%-8i                                  ->NVEGSEGMENT  \n',in.nvegsegment);
   else %% if in.nvegsegment is not defined, use default value of 1. However, do not enforce nvegsegment to be 1.
       fprintf(fid,'%-8i                                  ->NVEGSEGMENT  \n',1);
   end

   if isfield(in, 'ibreaking')
       fprintf(fid,'%-8i                                  ->IBREAKING  \n',in.ibreaking);
   else %% if in.ibreaking is not defined, use default value of 1. However, do not enforce ibreaking to be 1.
       fprintf(fid,'%-8i                                  ->IBREAKING  \n',1);
   end
   fprintf(fid,'%-8i                                  ->IDISS  \n',in.idiss);
   fprintf(fid,'%-8i                                  ->IFv  \n',in.iFv);
end
fprintf(fid,'%11.4f                                ->DXC\n',in.dx);
fprintf(fid,'%11.4f                                ->GAMMA \n',in.gamma);
% fprintf(fid,'%-8i                                -f>IWEIBULL \n',in.iweibull);
% if in.iweibull == 1;
%    fprintf(fid,'%11.4f                                ->DIKETOE \n',in.diketoe);
% end

if in.iprofl==1;
  fprintf(fid,'%11.4f%11.4f%11.4f         ->D50 WF SG\n',[in.d50, in.wf, in.sg]);
  fprintf(fid,'%11.4f%11.4f%11.4f%11.4f              ->EFFB EFFF SLP\n',...
          [in.effb, in.efff, in.slp, in.slpot]);
  fprintf(fid,'%11.4f%11.4f                    ->TANPHI BLP\n',[in.tanphi, in.blp]);
end
if in.iover;
  fprintf(fid,'%11.4f                               ->RWH \n',in.rwh);
end
if in.iperm;
  fprintf(fid, '%11.4f%11.4f%11.4f\n',in.stoneporo, in.stonedia, in.criticalstability )
end

fprintf(fid,'%-8i                                  ->ILAB\n',in.ilab);
% if in.iprofl==0;
%   fprintf(fid,'%11.4f%11.4f%11.4f%11.4f\n',[in.Tp(i), in.Hrms(i), in.Wsetup(i), in.angle(i)]);
% else
if in.ilab==1
  fprintf(fid,'%-8i                                  ->NWAVE \n',in.nwave);
  fprintf(fid,'%-8i                                  ->NSURGE \n',in.nsurg);
  for i = 1:length(in.Hrms)
    if in.iveg<=2
            fprintf(fid,'%11.2f%11.4f%11.4f%11.4f%11.4f%11.4f\n', ...
            [in.timebc_wave(i) in.Tp(i), in.Hrms(i), in.Wsetup(i), in.swlbc(i), in.angle(i)]);
    elseif in.iveg==3
       if in.idiss==1
          fprintf(fid,'%11.2f%11.4f%11.4f%11.4f%11.4f%11.4f\n', ...
             [in.timebc_wave(i) in.Tp(i), in.Hrms(i), in.Wsetup(i), in.swlbc(i), in.angle(i)]);
       elseif in.idiss>=2
          fprintf(fid,'%11.2f%11.4f%11.4f%11.4f%11.4f%11.4f%11.4f%11.4f%11i%11.2f\n', ...
            [in.timebc_wave(i), in.Tp(i), in.Hrms(i), in.Wsetup(i), in.swlbc(i), in.angle(i),...
             in.freqmin(i),in.freqmax(i),in.numfreq(i),in.JONSWAPgamma(i)]);
       end
    end
  end
else
  fprintf(fid,'%-8i                                  ->NWAVE \n',in.nwave-1);
  fprintf(fid,'%-8i                                  ->NSURGE \n',in.nsurg-1);
  for i = 1:length(in.Hrms)
    if in.iveg<=2
      fprintf(fid,'%11.2f%11.4f%11.4f%11.4f\n', ...
      [in.timebc_wave(i) in.Tp(i), in.Hrms(i),in.angle(i)]);
    elseif in.iveg==3
       if in.idiss==1
          fprintf(fid,'%11.2f%11.4f%11.4f%11.4f\n', ...
                [in.timebc_wave(i) in.Tp(i), in.Hrms(i),in.angle(i)]);
       elseif in.idiss==2
          fprintf(fid,'%11.2f%11.4f%11.4f%11.4f%11.4f%11.4f%11i%11.4f\n', ...
                [in.timebc_wave(i), in.Tp(i), in.Hrms(i),in.angle(i),...
                 in.freqmin(i), in.freqmax(i), in.numfreq(i),in.JONSWAPgamma(i)]);
        end
    end 
  end
  
  for i = 1:length(in.swlbc)
    fprintf(fid,'%11.2f%11.4f\n', ...
            [in.timebc_surg(i) in.swlbc(i)]);
  end

end
fprintf(fid,'%-8i                             ->NBINP \n',length(in.x));


% 
if ~isfield(in, 'x_p')
    in.x_p = in.x ; 
    in.zb_p = in.zb ; 
end

if in.iperm==1|| in.isedav>=1
  fprintf(fid,'%-8i                             ->NPINP \n',length(in.x_p));
end

dum = [in.x(:) in.zb(:) in.fw(:)];
fprintf(fid,'%11.6f %11.6f %11.6f\n',dum');

if in.iperm==1|in.isedav>=1
  dum = [in.x_p(:) in.zb_p(:)];
  fprintf(fid,'%11.4f%11.4f\n',dum');
end

%% prepare veg inputs
ind = find(in.x>=max(in.x)*in.veg_extent(1, 1)&in.x<=max(in.x)*in.veg_extent(1, 2));

if in.iveg==1
    dum = in.veg_Cd' ; 
    fprintf(fid,'%11.11f \n',dum');
    
    dum = [in.veg_n', in.veg_dia', in.veg_ht', in.veg_rod'];
    fprintf(fid,'%11.6f %11.6f %11.6f %11.6f\n',dum');   

elseif in.iveg==3 
    % only prepare veg Cd and Cdm when ivegCd = 0
    if in.ivegCd==0
        dum = [in.veg_Cd', in.veg_Cdm'] ; 
        fprintf(fid,'%11.11f %11.11f  \n',dum');
    end

    if in.ivegtype==0
        dum = [in.veg_n', in.veg_nblade', in.veg_dia', in.veg_diablade', in.veg_ht', in.veg_htblade', in.veg_thicknessblade', in.veg_rod'];
        fprintf(fid,'%11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f\n',dum');
    elseif in.ivegtype ==1       
        dum = [in.veg_n', in.veg_nblade', in.veg_dia', in.veg_diablade', in.veg_ht', in.veg_htblade', in.veg_thicknessblade', in.veg_Estem', in.veg_Eblade', in.veg_rod'];
        fprintf(fid,'%11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %11.6f %15.10f %15.10f %11.6f \n',dum');
        
    end
    
end

if in.iwind==1 
    if in.ilab==1
        fprintf(fid,'%-8i                                  ->NWIND \n',in.nwind);
        for i = 1:length(in.wind10)   
            fprintf(fid,'%11.4f%11.4f%11.4f\n', [in.time_wind(i) in.wind10(i), in.windangle(i)]);
        end
    end
end

fclose(fid);

%disp('Finished writing the infile')
