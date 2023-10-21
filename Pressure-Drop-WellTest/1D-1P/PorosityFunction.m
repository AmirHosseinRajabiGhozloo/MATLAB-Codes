function Phi = PorosityFunction(P,N,P0,Phi0,Cphi);

phi = @(p)[Phi0*(1+Cphi*(p-P0))];

for i = 1:N
    
    Phi(i) = phi(P(i));
    
end
end