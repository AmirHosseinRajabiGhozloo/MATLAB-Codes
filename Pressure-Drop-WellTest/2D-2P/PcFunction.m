function Pc = PcFunction(sWr,sOr,ndx,ndy,sw)

swn = zeros(ndx,ndy);
Pc = zeros(ndx,ndy,2);

for i = 1:ndx
    for j = 1:ndy 

        swn(i,j)=(sw(i,j)-sWr)/(1-sWr-sOr);

        if swn(i,j) <= 0.001

            Pc(i,j,1) = 150;

            Pc(i,j,2) = 15*(-1/3)*(1/(1-sOr-sWr))/0.001^(4/3);

        else

            Pc(i,j,1) = 15/swn(i,j)^(1/3);

            Pc(i,j,2) = 15*(-1/3)*(1/(1-sOr-sWr))/swn(i,j)^(4/3);

        end
    end
end
end