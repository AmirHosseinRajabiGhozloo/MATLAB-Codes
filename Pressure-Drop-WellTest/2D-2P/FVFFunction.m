function B = FVFFunction(Po,Pw,ndx,ndy,Bo0,Bw0,co,cw,P0)

FVF_oil = @(t)[Bo0/(1+co*(t-P0))];
FVF_wat = @(t)[Bw0/(1+cw*(t-P0))];

B = zeros(ndx,ndy,4);

for i = 1:ndx
    for j = 1:ndy

        B(i,j,1) = FVF_oil(Po(i,j));
        B(i,j,2) = FVF_wat(Pw(i,j));
        B(i,j,3) = co/Bo0;
        B(i,j,4) = cw/Bw0;

    end
end
end