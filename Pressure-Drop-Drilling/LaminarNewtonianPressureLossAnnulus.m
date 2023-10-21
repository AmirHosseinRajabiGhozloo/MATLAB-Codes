function [pl]=LaminarNewtonianPressureLossAnnulus(v,do,din,teta300)
mu=teta300;
pl=(mu*v)/(1000*((do-din)^2));
end