function [pl]=BinghamTurbulentPressureLossAnnulus(rho,v,teta300,teta600,do,din)
mu=teta600-teta300;
de=0.816*(do-din);
pl=((rho^0.75)*(v^1.75)*(mu^0.25))/(1800*(de^1.25));
end