function w=vfall(d,T,S)
%VFALL is a function to estimate fall velocity based on Soulsby's (1997)
%      optimization.
%      w = kvis/d* [sqrt(10.36^2 + 1.049 D^3) - 10.36]
%
% SYNTAX:  w=vfall(d,T,S)
%  where w = sediment fall speed (m/s)
%        d = grain diameter (mm)
%        T = temperature (deg C)
%        S = Salinity (ppt)
%
% Jarrell Smith
% Coastal and Hydraulics Laboratory
% Engineer Research and Development Center
% Vicksburg, MS

g=9.81;
rho=denfun(T,S);
kvis=kvisfun(T);
rhos=2650;
d=d/1000; %convert mm to m
s=rhos/rho;
D=(g*(s-1)/kvis^2)^(1/3)*d;
w=kvis./d.*(sqrt(10.36^2+1.049*D.^3)-10.36); %settling speed
