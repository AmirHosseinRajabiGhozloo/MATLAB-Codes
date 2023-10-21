function [Re]=NewtonianReynoldsAnnulus(rho,v,din,do,teta300)
mu=teta300;
de=0.816*(do-din);
Re=(928*rho*v*de)/(mu);
end