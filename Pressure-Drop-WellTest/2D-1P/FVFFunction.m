function Bo = FVFFunction(P,Nx,Ny,Co,B0,P0,BCN,BCE,BCS,BCW);

FVF = @(p)[B0/(1+Co*(p-P0))];

for j = 1:Ny
    
    for i = 1:Nx
    
        if BCN(i,j) == 100
       
            Bo(i,j,1) = FVF((P(i,j)+P(i-1,j))/2);
    
        else
            
            Bo(i,j,1) = FVF(P(i,j));
   
        end
        
        if BCE(i,j) == 100
        
            Bo(i,j,2) = FVF((P(i,j)+P(i,j+1))/2);
   
        else
            
            Bo(i,j,2) = FVF(P(i,j));
   
        end
        
        if BCS(i,j) == 100
       
            Bo(i,j,3) = FVF((P(i,j)+P(i+1,j))/2);
   
        else
            
            Bo(i,j,3) = FVF(P(i,j));
    
        end
        
        if BCW(i,j) == 100
       
            Bo(i,j,4) = FVF((P(i,j)+P(i,j-1))/2);
   
        else
            
            Bo(i,j,4) = FVF(P(i,j));
   
        end
        
        Bo(i,j,5) = FVF(P(i,j));
   
    end
end
end










