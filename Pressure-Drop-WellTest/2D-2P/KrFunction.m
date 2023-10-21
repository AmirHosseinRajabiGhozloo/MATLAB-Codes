function kr = KrFunction(sWr,sOr,ndx,ndy,sw)

swn = zeros(ndx,ndy);
kr = zeros(ndx,ndy,2);

for i = 1:ndx
    for j = 1:ndy

        swn(i,j) = (sw(i,j)-sWr)/(1-sWr-sOr);
        kr(i,j,1) = 0.9*(1-swn(i,j))^2;
        kr(i,j,2) = 0.6*swn(i,j)^2;

    end
end
end