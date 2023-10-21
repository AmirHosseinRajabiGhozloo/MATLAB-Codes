function phi = PorosityFunction(P,cphi,phi0,P0,ndx,ndy)

phi1 = @(t,l)[l*(1+cphi*(t-P0))];
phi = zeros(ndx,ndy);

for i = 1:ndx
    for j = 1:ndy

        phi(i,j) = phi1(P(i,j),phi0(i,j)); 

    end
end
end