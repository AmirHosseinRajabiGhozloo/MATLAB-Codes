function miu = ViscosityFunction(P,Nx,Ny,BCN,BCE,BCS,BCW);

viscosity = @(p)[2+2.5*10^(-5)*(p-14.7)];

for j = 1:Ny
    
    for i = 1:Nx
    
        if BCN(i,j) == 100
        
            miu(i,j,1) = viscosity((P(i,j)+P(i-1,j))/2);
    
        else
            
            miu(i,j,1) = viscosity(P(i,j));
    
        end
        
        if BCE(i,j) == 100
       
            miu(i,j,2) = viscosity((P(i,j)+P(i,j+1))/2);
    
        else
            
            miu(i,j,2) = viscosity(P(i,j));
    
        end
        
        if BCS(i,j) == 100
       
            miu(i,j,3) = viscosity((P(i,j)+P(i+1,j))/2);
    
        else
            
            miu(i,j,3) = viscosity(P(i,j));
    
        end
        
        if BCW(i,j) == 100
       
            miu(i,j,4) = viscosity((P(i,j)+P(i,j-1))/2);
    
        else
            
            miu(i,j,4) = viscosity(P(i,j));
    
        end
        
        miu(i,j,5) = viscosity(P(i,j));
    
    end
end
end










