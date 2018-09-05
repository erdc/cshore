function makeinfile_usace_vegfeature(in)
% function makeinfile_usace(in)

%First check for NaN
if max(isnan([in.x(:);in.zb(:);in.Tp(:);in.Hrms(:);in.swlbc(:);in.Wsetup(:);in.angle(:)]))
  error('The in structure contains NaN')
end


fid = fopen('infile','w');
fprintf(fid,'%i \n',length(in.header));
for i = 1:length(in.header)
  fprintf(fid,'%s \n',cell2mat(in.header(i)));
end
fprintf(fid,'%-8i                                  ->ILINE\n',in.iline);
fprintf(fid,'%-8i                                  ->IPROFL\n',in.iprofl);
if in.iprofl==1
  fprintf(fid,'%-8i                                  ->ISEDAV\n',in.isedav);
end
fprintf(fid,'%-8i                                  ->IPERM\n',in.iperm);
fprintf(fid,'%-8i                                  ->IOVER\n',in.iover);
if in.iover
  fprintf(fid,'%-8i                                  ->IWTRAN\n',in.iwtran);
  if in.iwtran==0
    fprintf(fid,'%-8i                                  ->IPOND\n',in.ipond);
  end
end
if in.iover==1&in.iperm==0&in.iprofl==1
  fprintf(fid,'%-8i                                  ->INFILT\n',in.infilt);
end
fprintf(fid,'%-8i                                  ->IWCINT\n',in.iwcint);
fprintf(fid,'%-8i                                  ->IROLL \n',in.iroll);
fprintf(fid,'%-8i                                  ->IWIND \n',in.iwind);
fprintf(fid,'%-8i                                  ->ITIDE \n',in.itide);
fprintf(fid,'%-8i                                  ->IVEG  \n',in.iveg);
fprintf(fid,'%-8i                                  ->IDISS  \n',in.idiss);
fprintf(fid,'%-8i                                  ->IFv  \n',in.iFv);
fprintf(fid,'%11.4f                                ->DXC\n',in.dx);
fprintf(fid,'%11.4f                                ->GAMMA \n',in.gamma);
fprintf(fid,'%-8i                                ->IWEIBULL \n',in.iweibull);
if in.iweibull == 1;
   fprintf(fid,'%11.4f                                ->DIKETOE \n',in.diketoe);
end

if in.iprofl==1;
  fprintf(fid,'%11.4f%11.4f%11.4f         ->D50 WF SG\n',[in.d50, in.wf, in.sg]);
  fprintf(fid,'%11.4f%11.4f%11.4f%11.4f              ->EFFB EFFF SLP\n',...
          [in.effb, in.efff, in.slp, in.slpot]);
  fprintf(fid,'%11.4f%11.4f                    ->TANPHI BLP\n',[in.tanphi, in.blp]);
end
if in.iover;
  fprintf(fid,'%11.4f                               ->RWH \n',in.rwh);
end
fprintf(fid,'%-8i                                  ->ILAB\n',in.ilab);
% if in.iprofl==0;
%   fprintf(fid,'%11.4f%11.4f%11.4f%11.4f\n',[in.Tp(i), in.Hrms(i), in.Wsetup(i), in.angle(i)]);
% else
if in.ilab==1
  fprintf(fid,'%-8i                                  ->NWAVE \n',in.nwave);
  fprintf(fid,'%-8i                                  ->NSURGE \n',in.nsurg);
  for i = 1:length(in.Hrms)
    if in.idiss==1
        fprintf(fid,'%11.2f%11.4f%11.4f%11.4f%11.4f%11.4f\n', ...
            [in.timebc_wave(i) in.Tp(i), in.Hrms(i), in.Wsetup(i), in.swlbc(i), in.angle(i)]);
    elseif in.idiss>=2
        fprintf(fid,'%11.2f%11.4f%11.4f%11.4f%11.4f%11.4f%11.4f%11.4f%11i\n', ...
            [in.timebc_wave(i), in.Tp(i), in.Hrms(i), in.Wsetup(i), in.swlbc(i), in.angle(i),...
             in.freqmin(i),in.freqmax(i),in.numfreq(i)]);
    end
  end
else
  fprintf(fid,'%-8i                                  ->NWAVE \n',in.nwave-1);
  fprintf(fid,'%-8i                                  ->NSURGE \n',in.nsurg-1);
  for i = 1:length(in.Hrms)
    if in.idiss==1
        fprintf(fid,'%11.2f%11.4f%11.4f%11.4f\n', ...
                [in.timebc_wave(i) in.Tp(i), in.Hrms(i),in.angle(i)]);
    elseif in.idiss==2
        fprintf(fid,'%11.2f%11.4f%11.4f%11.4f%11.4f%11.4f%11.4f\n', ...
                [in.timebc_wave(i), in.Tp(i), in.Hrms(i),in.angle(i),...
                 in.freqmin(i), in.freqmax(i), in.numfreq(i)]);
    end
        
  end
  
  for i = 1:length(in.swlbc)
    fprintf(fid,'%11.2f%11.4f\n', ...
            [in.timebc_surg(i) in.swlbc(i)]);
  end

end
fprintf(fid,'%-8i                             ->NBINP \n',length(in.x));
dum = [in.x(:) in.zb(:) in.fw(:)];
fprintf(fid,'%11.6f%11.6f%11.6f\n',dum');
if in.iveg==1
    fprintf(fid,'%10.6f                                ->VEGCD\n',in.veg_Cd );
    fprintf(fid,'%10.6f                                ->VEGCDM\n',in.veg_Cdm );
    dum = zeros(length(in.x(:)),4);
    ind = find(in.x>=max(in.x)*in.veg_extent(1)&in.x<=max(in.x)*in.veg_extent(2));
    dum(ind,:) = repmat([in.veg_n in.veg_dia in.veg_ht in.veg_rod],length(ind),1);
    fprintf(fid,'%11.6f%11.6f%11.6f%11.6f\n',dum');
end

fclose(fid);

%disp('Finished writing the infile')
