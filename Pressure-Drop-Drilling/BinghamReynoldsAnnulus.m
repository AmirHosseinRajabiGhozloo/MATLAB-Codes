function [Re]=BinghamReynoldsAnnulus(rho,v,do,din,teta600,teta300)
mu=teta600-teta300;
ty=teta300-mu;
mua=mu+(5*ty*(do-din))/(v);
de=0.816*(do-din);
Re=(928*rho*v*de)/(mua);
end