function rho=denfun(T,S)
%DENFUN estimates water density from temperature and salinity
%Approximation from VanRijn, L.C. (1993) Handbook for Sediment Transport 
%                     by Currents and Waves
%
%SYNTAX:  rho=denfun(T,S)
%  where  rho = density of water (kg/m^3)
%           T = temperature (C)
%           S = salinity (o/oo)
%

T=reshape(T,length(T),1);
S=reshape(S,length(S),1);
CL=(S-0.03)/1.805;  %VanRijn
in=CL<0;
CL(in)=0;
rho=1000 + 1.455.*CL - 6.5e-3* (T-4+0.4.*CL).^2;  %from VanRijn (1993)
