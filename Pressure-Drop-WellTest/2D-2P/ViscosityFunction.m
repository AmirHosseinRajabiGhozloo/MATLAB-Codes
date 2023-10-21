function miu = ViscosityFunction(po,pw,ndx,ndy)

viscosity_oil = @(t)[1];                         %% cp
viscosity_wat = @(t)[1];                         %% cp

miu = zeros(ndx,ndy,2);

for i = 1:ndx
    for j = 1:ndy

        miu(i,j,1) = viscosity_oil(po(i,j));
        miu(i,j,2) = viscosity_wat(pw(i,j));

    end
end