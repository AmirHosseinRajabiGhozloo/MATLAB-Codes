function [Re]=NewtonianReynoldsPipe(rho,v,d,teta300)
mu=teta300;
Re=(928*rho*v*d)/(mu);
end