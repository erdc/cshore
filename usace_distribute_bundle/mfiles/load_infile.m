function in = load_infile(fname)
if nargin==0
    fname='infile';
end

disp(['loading info from file named ',fname])

fid = fopen([fname]);
tot = textscan(fid,'%s','delimiter','\n');
tot = tot{:};
fclose(fid);
%
nhlines = str2num(tot{1});
in.tot = tot;
in.header = {tot{2:1+nhlines}}';
%
keys = {'ILINE';'IPROFL';'ISEDAV';'IPERM';'IOVER';'IWTRAN';'IPOND';'INFILT';'IWCINT';
        'IROLL ';'IWIND ';'ITIDE ';'IVEG  ';'DXC';'GAMMA';'RWH';'ILAB';'NWAVE';'NSURGE';'NBINP'} ;
for i = 1:length(keys)
  row_ind = find(contains(tot,keys{i},'IgnoreCase',true));
  if ~isempty(row_ind)
    dum = split(convertCharsToStrings(tot{row_ind}));
    in=setfield(in,lower(keys{i}),str2num(dum{1}));
  end
end
in.dx=in.dxc;
in.nsurg=in.nsurge;
if in.iprofl>0
    row_ind = find(contains(tot,'D50','IgnoreCase',true));
    dum = split(convertCharsToStrings(tot{row_ind}));
    in.d50 = str2num(dum{1});
    in.wf = str2num(dum{2});
    in.sg = str2num(dum{3});
    row_ind = find(contains(tot,'EFFB','IgnoreCase',true));
    dum = split(convertCharsToStrings(tot{row_ind}));
    in.effb = str2num(dum{1});
    in.efff = str2num(dum{2});
    in.slp = str2num(dum{3});
    in.slpot = str2num(dum{4});
    row_ind = find(contains(tot,'TANPHI','IgnoreCase',true));
    dum = split(convertCharsToStrings(tot{row_ind}));
    in.tanphi = str2num(dum{1});
    in.blp = str2num(dum{2});
end

if in.ilab==0
    row_ind = find(contains(tot,'NWAVE','IgnoreCase',true))+2;
    dum = cell2mat(cellfun(@str2num,({tot{row_ind:row_ind+in.nwave}}'),'UniformOutput',false));
    in.timebc_wave = dum(:,1)';
    in.Tp= dum(:,2)';
    in.Hrms= dum(:,3)';
    in.Wsetup= zeros(size(in.Hrms));
    in.angle = dum(:,4)';
    dum = cell2mat(cellfun(@str2num,({tot{row_ind+in.nwave+1:row_ind+in.nwave+in.nsurge+1}}'),'UniformOutput',false));
    in.timebc_surg = dum(:,1)';
    in.swlbc= dum(:,2)';
    in.nwave = in.nwave+1;
    in.nsurg = in.nsurg+1;
elseif in.ilab==1
    row_ind = find(contains(tot,'NWAVE','IgnoreCase',true))+2;
    dum = cell2mat(cellfun(@str2num,({tot{row_ind:row_ind+in.nwave}}'),'UniformOutput',false));
    in.timebc_wave = dum(:,1)';
    in.Tp= dum(:,2)';
    in.Hrms= dum(:,3)';
    in.Wsetup= dum(:,4)';
    in.swlbc= dum(:,5)';
    in.angle = dum(:,6)';
end

row_ind = find(contains(tot,'NBINP','IgnoreCase',true))+1;
dum = cell2mat(cellfun(@str2num,({tot{row_ind:row_ind+in.nbinp-1}}'),'UniformOutput',false));
in.x= dum(:,1)';
in.zb= dum(:,2)';
in.fw= dum(:,3)';





return

% find ILINE 
row_ind = find(contains(tot,'iline','IgnoreCase',true))
% find IPROFL 
row_ind = find(contains(tot,'iprofl','IgnoreCase',true))
dum = split(convertCharsToStrings(tot{row_ind}));
in.iprofl = str2num(dum{1});



