function kvis=kvisfun(T)
%KVISFUN estimates kinematic viscosity of water
%Approximation from VanRijn, L.C. (1989) Handbook of Sediment Transport
%
%SYNTAX:  kvis=kvisfun(T)
%  where  kvis = kinematic viscosity (m^2/sec)
%            T = temperature (C)
%

kvis=1e-6*(1.14 - 0.031*(T-15) + 6.8E-4*(T-15).^2);
