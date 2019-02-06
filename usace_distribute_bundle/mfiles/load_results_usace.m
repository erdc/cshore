function out = load_results_usace(flg)
% function out = load_results_usace(flg)

fid = fopen('ODOC');
tot = textscan(fid,'%s','delimiter','\n');
tot = tot{:};
fclose(fid);

% find header
dum =strfind(tot,'OPTION ILINE');
row_ind = find(~cellfun('isempty',dum));
if isempty(row_ind)
  out.run_success = 0;
  out.header = tot;
  return
end
out.header = tot(1:row_ind-1);

% find IPROFL 
dum =strfind(tot,'OPTION IPROFL');
row_ind = find(~cellfun('isempty',dum));
row = tot{row_ind};
col_ind = find(row=='=');
out.params.iprofl=str2num(row(col_ind+1:end));

% find ISEDAV 
dum =strfind(tot,'ISEDAV');
row_ind = find(~cellfun('isempty',dum));
if ~isempty(row_ind)
  row = tot{row_ind};
  col_ind = find(row=='=');
  out.params.isedav=str2num(row(col_ind+1:col_ind+3));
else
  out.params.isedav=0;
end



% find IPERM
dum =strfind(tot,'IMPERMEABLE');
row_ind = find(~cellfun('isempty',dum));
out.params.iperm=isempty(row_ind);

% find NBINP
dum =strfind(tot,'NBINP');
row_ind = find(~cellfun('isempty',dum));
row=cell2mat(tot(row_ind));
col_ind = find(row=='=');
out.params.nbinp=str2num(row(col_ind(end)+1:end));

% find GAMMA (all GAMMA)
dum =strfind(tot,'Gamma')';
row_ind = find(~cellfun('isempty',dum));
row = tot(row_ind);
col_ind = find(cell2mat(row(1))=='=');
row = cell2mat(row);
gamma=str2num(row(:,col_ind+1:end));

% get longshore transport
dum =strfind(tot,'Transport Rate');
row_ind = find(~cellfun('isempty',dum));
for k = 1:length(row_ind)
  row = tot{row_ind(k)};
  col_ind = find(row=='=');
  out.sed.longshore_transport(k)=str2num(row(col_ind(end)+1:end));
end

%get wave conditions at SB
dum =strfind(tot,'INPUT WAVE');
row_ind = find(~cellfun('isempty',dum));
ind_begin =  row_ind+3;
dum =strfind(tot,'INPUT BEACH AND STRUCTURE');
row_ind = find(~cellfun('isempty',dum));
ind_end = row_ind-2;
cnt = 0;
wave_cond = [];
for i = ind_begin:ind_end
  cnt = cnt+1;
  str2num(cell2mat(tot(i,:)));
  wave_cond=[wave_cond;str2num(cell2mat(tot(i,:)))];
end

out.bc.time_offshore = wave_cond(:,1);
out.bc.Tp_offshore=wave_cond(:,2);
out.bc.Hrms_offshore=wave_cond(:,3);
out.bc.wave_setup_offshore=wave_cond(:,4);
out.bc.strm_tide_offshore=wave_cond(:,5);
out.bc.angle_offshore=wave_cond(:,6);


% find  runup
dum2p   = strfind(tot,'2 percent runup');row_ind = ~cellfun('isempty',dum2p);  dum2p = cell2mat(tot(row_ind));
dummean = strfind(tot,'Mean runup');     row_ind = ~cellfun('isempty',dummean);dummean = cell2mat(tot(row_ind));
if ~isempty(dum2p)
  for ii = 1:size(dum2p,1)
    col_ind = strfind(dum2p(ii,:),'R2P=');
    dum22p =     str2num(dum2p(ii,col_ind+4:end));
    dum2mean = str2num(dummean(ii,col_ind+4:end));
    if ~isempty(dum22p)
      runup_2_percent(ii)=dum22p;
      runup_mean(ii)     =dum2mean;
    else
      runup_2_percent(ii)=NaN;
      runup_mean(ii)=NaN;
    end
    out.hydro(ii).runup_2_percent = runup_2_percent(ii);
    out.hydro(ii).runup_mean = runup_mean(ii);
  end
end

% find  jdry
dum   = strfind(tot,'JDRY');row_ind = ~cellfun('isempty',dum);  dum = cell2mat(tot(row_ind));
if ~isempty(dum)
  for ii = 1:size(dum,1)
    col_ind = strfind(dum(ii,:),'JDRY=');
    dum2 =    str2num(dum(ii,col_ind+5:end));
    if ~isempty(dum2)
      jdry(ii)=dum2;
    else
      jdry(ii)=NaN;
    end
    out.hydro(ii).jdry = jdry(ii);
  end
end

% find  SWL at sea boundary
dum   = strfind(tot,' SWL=');row_ind = ~cellfun('isempty',dum);  dum = cell2mat(tot(row_ind));
if ~isempty(dum)
  for ii = 1:size(dum,1)
    col_ind = strfind(dum(ii,:),'SWL=');
    dum2 =    str2num(dum(ii,col_ind+5:end));
    if ~isempty(dum2)
      swl(ii)=dum2;
    else
      swl(ii)=NaN;
    end
    out.hydro(ii).swl = swl(ii);
  end
end

% find  node number of SWL 
dum   = strfind(tot,' JSWL=');row_ind = ~cellfun('isempty',dum);  dum = cell2mat(tot(row_ind));
if ~isempty(dum)
  for ii = 1:size(dum,1)
    col_ind = strfind(dum(ii,:),'JSWL=');
    dum2 =    str2num(dum(ii,col_ind+5:end));
    if ~isempty(dum2)
      jswl(ii)=dum2;
    else
      jswl(ii)=NaN;
    end
    out.hydro(ii).jswl = jswl(ii);
  end
end

% find  jr
dum   = strfind(tot,'JR=');row_ind = ~cellfun('isempty',dum);dum = cell2mat(tot(row_ind));
if ~isempty(dum)
  for ii = 1:size(dum,1)
    col_ind = strfind(dum(ii,:),'JR=');
    dum2 =    str2num(dum(ii,col_ind+5:end));
    if ~isempty(dum2)
      jr(ii)=dum2;
    else
      jr(ii)=NaN;
    end
    out.hydro(ii).jr = jr(ii);
  end
end

% swash zone bottom slope
dum =strfind(tot,'Swash zone bottom slope');
row_ind = ~cellfun('isempty',dum);
row_ind = find(row_ind);
dum_slp = cell2mat(tot(row_ind));
dum_x1  = cell2mat(tot(row_ind+1));
dum_x2  = cell2mat(tot(row_ind+2));
dum_z1  = cell2mat(tot(row_ind+3));
dum_z2  = cell2mat(tot(row_ind+4));
if ~isempty(dum_slp)
  for ii = 1:size(dum_slp,1)
    col_ind = strfind(dum_slp(ii,:),'=');
    dum_slp2 = str2num(dum_slp(ii,col_ind+1:end));
    dum_x12 = str2num(dum_x1(ii,col_ind+1:end));
    dum_x22 = str2num(dum_x2(ii,col_ind+1:end));
    dum_z12 = str2num(dum_z1(ii,col_ind+1:end));
    dum_z22 = str2num(dum_z2(ii,col_ind+1:end));
    if ~isempty(dum_slp2)
      slprun(ii)=dum_slp2;
      x1run(ii)=dum_x12;
      x2run(ii)=dum_x22;
      z1run(ii)=dum_z12;
      z2run(ii)=dum_z22;
    else
      slprun(ii)=NaN;
      x1run(ii)=NaN;
      x2run(ii)=NaN;
      z1run(ii)=NaN;
      z2run(ii)=NaN;
    end
    out.hydro(ii).slprun = slprun(ii);    
    out.hydro(ii).x1run = x1run(ii);
    out.hydro(ii).x2run = x2run(ii);
    out.hydro(ii).z1run = z1run(ii);
    out.hydro(ii).z2run = z2run(ii);
  end
end

% %%%%%%%%%%%%Get info from the infile%%%%%%%%%%%%%%%%%%
fid = fopen('infile');
tot = textscan(fid,'%s','delimiter','\n');
tot = tot{:};

% find iover 
dum =strfind(tot,'IOVER');
row_ind = ~cellfun('isempty',dum);
dum = cell2mat(tot(row_ind));
col_ind = strfind(dum,'-');
out.params.iover = str2num(dum(1:col_ind-1));

% find iveg 
dum =strfind(tot,'IVEG');
row_ind = ~cellfun('isempty',dum);
dum = cell2mat(tot(row_ind));
col_ind = strfind(dum,'-');
out.params.iveg = str2num(dum(1:col_ind-1));

% find effB and effF and blp
if out.params.iprofl
  dum =strfind(tot,'EFFB');
  row_ind = ~cellfun('isempty',dum);
  dum = cell2mat(tot(row_ind));
  if ~isempty(dum)
    col_ind = strfind(dum,'-');
    dum = str2num(dum(1:col_ind-1));
    out.params.effB=dum(1);
    out.params.effF=dum(2);
  end

  % find blp and tanphi
  dum =strfind(tot,'BLP');
  row_ind = ~cellfun('isempty',dum);
  dum = cell2mat(tot(row_ind));
  col_ind = strfind(dum,'-');
  dum = str2num(dum(1:col_ind-1));
  out.params.tanphi=dum(1);
  out.params.blp=dum(2);
  % find ilab
  dum =strfind(tot,'ILAB');
  row_ind = ~cellfun('isempty',dum);
  dum = cell2mat(tot(row_ind));
  col_ind = strfind(dum,'-');
  out.params.ilab = str2num(dum(1:col_ind-1));
end
% find vegitation extent
if out.params.iveg
  dum =strfind(tot,'VEGCD');
  row_ind = find(~cellfun('isempty',dum));
  dum = tot(row_ind+1:row_ind+out.params.nbinp);
  dum=cell2mat(cellfun(@str2num,dum,'UniformOutput',false));
  out.veg.n=dum(:,1);
  out.veg.dia=dum(:,2);
  out.veg.ht=dum(:,3);
  out.veg.rod=dum(:,4);
end



fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid=fopen('OBPROF');
cnt=0;
while 1
  tline = fgetl(fid);
  if ~ischar(tline), break, end
  cnt = cnt+1;
  tline = str2num(tline);
  if length(tline)==3;
    N = tline(2);
    tme = tline(3);
  elseif length(tline)==2;
    N = tline(1);
    tme=tline(2);
  end
  if (out.params.iveg&cnt>1)&out.params.isedav==0
    [tot]=fscanf(fid,'%f %f %f\n',[3,N])';
    out.morpho(cnt).ivegitated = tot(:,3);
  elseif out.params.isedav==1
    [tot]=fscanf(fid,'%f %f %f\n',[3,N])';
    out.morpho(cnt).zb_p = tot(:,3);
  else
    [tot]=fscanf(fid,'%f %f \n',[2,N])';
  end
  out.morpho(cnt).time = tme;
  out.morpho(cnt).x = tot(:,1);
  out.morpho(cnt).zb = tot(:,2);
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen('OSETUP');
%	  WRITE(22,1500) XB(J),(WSETUP(J)+SWLBC(IWAVE)),H(J),SIGMA(J)
cnt=0;
while 1
  tline = fgetl(fid);
  if ~ischar(tline), break, end
  cnt = cnt+1;
  tline = str2num(tline);
  if tline(1)==1
    N = tline(2);tme=tline(end);
  else
    N = tline(1);
  end
  [tot]=fscanf(fid,'%f %f %f %f \n',[4,N])';
  out.hydro(cnt).time_end = tme;
  out.hydro(cnt).x     = [tot(:,1); NaN(length(out.morpho(1).x)-size(tot,1),1)];
  out.hydro(cnt).setup = [tot(:,2); NaN(length(out.morpho(1).x)-size(tot,1),1)];
  out.hydro(cnt).depth = [tot(:,3); NaN(length(out.morpho(1).x)-size(tot,1),1)];
  out.hydro(cnt).sigma = [tot(:,4); NaN(length(out.morpho(1).x)-size(tot,1),1)];
  out.hydro(cnt).Hrms = sqrt(8)*out.hydro(cnt).sigma ;
end
fclose(fid);
num_output = cnt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen('OXVELO');
%    WRITE(27,1500) XB(J),UMEAN(J),USTD(J)
for i = 1:num_output
  tline = fgetl(fid);
  tline = str2num(tline);
  [tot]=fscanf(fid,'%f %f %f\n',[3,tline(2)])';
  out.hydro(i).x_xvelo = [tot(:,1); NaN(length(out.morpho(1).x)-size(tot,1),1)];
  out.hydro(i).umean = [tot(:,2); NaN(length(out.morpho(1).x)-size(tot,1),1)];
  out.hydro(i).ustd = [tot(:,3); NaN(length(out.morpho(1).x)-size(tot,1),1)];
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen('OYVELO');
%    WRITE(28,1500) XB(J),STHETA(J),VMEAN(J),VSTD(J)
for i = 1:num_output
  tline = fgetl(fid);
  if ~ischar(tline), break, end
  tline = str2num(tline);
  [tot]=fscanf(fid,'%f %f %f %f\n',[4,tline(2)])';
  if ~isempty(tot)&size(tot,1)>10
    out.hydro(i).x_yvelo = [tot(:,1); NaN(length(out.morpho(1).x)-size(tot,1),1)];
    out.hydro(i).stheta = [tot(:,2); NaN(length(out.morpho(1).x)-size(tot,1),1)];
    out.hydro(i).vmean = [tot(:,3); NaN(length(out.morpho(1).x)-size(tot,1),1)];
    out.hydro(i).vstd = [tot(:,4); NaN(length(out.morpho(1).x)-size(tot,1),1)];
  end

end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if out.params.iprofl
  fid=fopen('OCROSS');
  %    WRITE(32,1500) XB(J),QBX(J),QSX(J),(QBX(J)+QSX(J))
  for i = 1:num_output
    tline = fgetl(fid);
    tline = str2num(tline);
    [tot]=fscanf(fid,'%f %f %f %f\n',[4,tline(2)])';
    out.sed(i).x_cross = [tot(:,1); NaN(length(out.morpho(1).x)-size(tot,1),1)];
    out.sed(i).qbx = [tot(:,2); NaN(length(out.morpho(1).x)-size(tot,1),1)];
    out.sed(i).qsx = [tot(:,3); NaN(length(out.morpho(1).x)-size(tot,1),1)];
    out.sed(i).qx = [tot(:,4); NaN(length(out.morpho(1).x)-size(tot,1),1)];
  end
  fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if out.params.iprofl&0
  fid=fopen('OLONGS')
  %    WRITE(33,1500) XB(J),QBY(J),QSY(J),(QBY(J) + QSY(J))
  for i = 1:num_output
    tline = fgetl(fid)
    tline = str2num(tline);
    [tot]=fscanf(fid,'%f %f %f %f\n',[4,tline(2)])';
    out.sed(i).x_long = [tot(:,1); NaN(length(out.morpho(1).x)-size(tot,1),1)];
    out.sed(i).qby = [tot(:,2); NaN(length(out.morpho(1).x)-size(tot,1),1)];
    out.sed(i).qsy = [tot(:,3); NaN(length(out.morpho(1).x)-size(tot,1),1)];
    out.sed(i).qy = [tot(:,4); NaN(length(out.morpho(1).x)-size(tot,1),1)];
  end
  fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if out.params.iprofl
  fid=fopen('OBSUSL');
  %        IF(IPROFL.EQ.1) WRITE(30,1500) XB(J),PB(J),PS(J),VS(J)
  for i = 1:num_output
    tline = fgetl(fid);
    tline = str2num(tline);
    [tot]=fscanf(fid,'%f %f %f %f \n',[4,tline(2)])';
    out.sed(i).x_susl = [tot(:,1); NaN(length(out.morpho(1).x)-size(tot,1),1)];
    out.sed(i).ps = [tot(:,2); NaN(length(out.morpho(1).x)-size(tot,1),1)];
    out.sed(i).pb = [tot(:,3); NaN(length(out.morpho(1).x)-size(tot,1),1)];
    out.sed(i).vs = [tot(:,4); NaN(length(out.morpho(1).x)-size(tot,1),1)];
  end
  fclose(fid);
end

if isfield(out,'hydro');
  if (length(out.bc.Hrms_offshore)==length(out.hydro));  
    out.run_success =  1;
  else
    out.run_success = .5;
  end
else
  out.run_success = 0;
end

out = orderfields(out);

