function [pl]=LaminarNewtonianPressureLossPipe(v,d,teta300)
mu=teta300;
pl=(mu*v)/(1500*(d^2));
end
