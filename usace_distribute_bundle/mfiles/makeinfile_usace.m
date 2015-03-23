function makeinfile_2014(in)
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
fprintf(fid,'%11.4f                                ->DXC\n',in.dx);
fprintf(fid,'%11.4f                                ->GAMMA \n',in.gamma);
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
    fprintf(fid,'%11.2f%11.4f%11.4f%11.4f%11.4f%11.4f\n', ...
            [in.timebc_wave(i) in.Tp(i), in.Hrms(i), in.Wsetup(i), in.swlbc(i), in.angle(i)]);
  end
else
  fprintf(fid,'%-8i                                  ->NWAVE \n',in.nwave-1);
  fprintf(fid,'%-8i                                  ->NSURGE \n',in.nsurg-1);
  for i = 1:length(in.Hrms)
    fprintf(fid,'%11.2f%11.4f%11.4f%11.4f\n', ...
            [in.timebc_wave(i) in.Tp(i), in.Hrms(i),in.angle(i)]);
  end
  for i = 1:length(in.swlbc)
    fprintf(fid,'%11.2f%11.4f\n', ...
            [in.timebc_surg(i) in.swlbc(i)]);
  end

end
%end
fprintf(fid,'%-8i                             ->NBINP \n',length(in.x));
% if isfield(in,'xp')
% fprintf(fid,'%-8i                             ->NPINP \n',length(in.xp));
% else
%   fprintf(fid,'0                                    ->NPINP \n');
% end
dum = [in.x(:) in.zb(:) in.fw(:)];
fprintf(fid,'%11.4f%11.4f%11.4f\n',dum');
fclose(fid);