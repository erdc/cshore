boundary_gauge=1
load (['./waves_short'])  
load (['./profile'])  
offset = ft2m(10+3);

  
Hrms_data = [ 0.28 0.33 0.33  0.39  0.42];  
T_data = [4 4 6 4 6];
repeats = [6 18 18 15 64+23];
Hrms = [];
T = [];
for i = 1:length(Hrms_data)
  Hrms = [Hrms ;repmat(Hrms_data(i),repeats(i),1)];
  T = [T ;repmat(T_data(i),repeats(i),1)];
end
dt = 1/4;nwave=36/dt;dum = ones(1,nwave);
timebc = dt*3600*(cumsum(dum)-1);

fid = fopen('Hrms_bc.dat','w');
dum = [timebc(:) Hrms(:)];
fprintf(fid,'%10.2f %5.3f \n ',dum');
fclose(fid);
fid = fopen('T_bc.dat','w');
dum = [timebc(:) T(:)];
fprintf(fid,'%10.2f %5.3f \n ',dum');
fclose(fid);
%init bathy
x0 = waves_short(1).x_wg(boundary_gauge);
dx = 1;
zb =[profile(1).z]-offset;
x = [profile(1).x;];
x2 =x0:dx:max(x);
zb2 = interp1([0; x],[zb(1);zb],x2);
x2 = x2-x0;
fid = fopen('zb_init.dat','w');
dum = [x2(:) zb2(:)];
fprintf(fid,'%8.3f %8.3f \n ',dum');
fclose(fid);

fid = fopen('zb_time.dat','w');
fprintf(fid,'%8.3f \n ',[profile.time_under_wave_attack]*60);
fclose(fid);

fid = fopen('zb.dat','w');
zb = [profile.z]-offset;
fprintf(fid,[repmat('%8.3f',1,length(profile)),' \n '],zb');
fclose(fid);

fid = fopen('x.dat','w');
xall = [profile.x]-x0;
fprintf(fid,[repmat('%8.3f',1,length(profile)),' \n '],xall');
fclose(fid);


