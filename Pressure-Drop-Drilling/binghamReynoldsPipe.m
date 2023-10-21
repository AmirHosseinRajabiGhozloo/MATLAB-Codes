function [Re]=binghamReynoldsPipe(rho,v,d,teta600,teta300)
mu=teta600-teta300;
ty=teta300-mu;
mua=mu+(6.66*ty*d)/(v);
Re=(928*rho*v*d)/(mua);
end