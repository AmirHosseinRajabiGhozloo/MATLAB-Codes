function [pl]=LaminarBinghamPressureLossPipe(d,v,teta300,teta600)
mu=teta600-teta300;
ty=teta300-mu;
pl=(mu*v)/(1500*(d^2))+(ty)/(225*d);
end