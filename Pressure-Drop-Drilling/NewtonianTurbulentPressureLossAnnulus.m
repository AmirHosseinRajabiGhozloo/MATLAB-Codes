function [pl]=NewtonianTurbulentPressureLossAnnulus(rho,v,din,do,teta300)
mu=teta300;
de = 0.816 * (do-din);
pl=((rho^0.75)*(mu^0.25)*(v^1.75))/((1800)*(de^1.25));
end