function [pl]=BinghamTurbulentPressureLossPipe(rho,v,teta300,teta600,d)
mu=teta600-teta300;
pl=((rho^0.75)*(v^1.75)*(mu^0.25))/(1800*(d^1.25));
end