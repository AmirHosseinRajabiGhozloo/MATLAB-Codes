function [pl] = LaminarBinghamPressureLossAnnulus(din,do,v,teta300,teta600)
mu=teta600-teta300;
ty=teta300-mu;
pl=(mu*v)/(1000*((do-din)^2))+(ty)/(200*(do-din));
end
