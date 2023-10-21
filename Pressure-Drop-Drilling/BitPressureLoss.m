function [bpl]=BitPressureLoss(rho,q,cd,db1,db2,db3)
d1=db1/32;
d2=db2/32;
d3=db3/32;
At=(pi*((d1^2)+(d2^2)+(d3^2)))/4;
bpl=(8.311*10^(-5)*rho*q^2)/((cd^2)*(At^2));
end