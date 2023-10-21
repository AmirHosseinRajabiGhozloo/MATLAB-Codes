function [pl]=NewtonianTurbulentPressureLossPipe(rho,v,d,teta300)
mu=teta300;
pl=((rho^0.75)*(mu^0.25)*(v^1.75))/((1800)*(d^1.25));
end