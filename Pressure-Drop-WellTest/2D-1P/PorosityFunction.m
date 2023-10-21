function Phi = PorosityFunction(P,Nx,Ny,P0,Cphi,Phi0,BCN,BCE,BCS,BCW);

phi = @(p,Phi0)[Phi0*(1+Cphi*(p-P0))];

for j = 1:Ny
    
    for i = 1:Nx
           
        if BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block
        
            Phi(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block
        
            Phi(i,j) = 0;
            
        elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block
        
            Phi(i,j) = 0;
            
        elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
        
            Phi(i,j) = 0;
            
        elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block
        
            Phi(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block
        
            Phi(i,j) = 0;
            
        elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block
        
            Phi(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block
        
            Phi(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block
        
            Phi(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
        
            Phi(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block
        
            Phi(i,j) = 0;
            
        elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block
        
            Phi(i,j) = 0;
            
        elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block
        
            Phi(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
        
            Phi(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block
        
            Phi(i,j) = 0;
            
        elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
        
            Phi(i,j) = 0;
            
        else
            
            Phi(i,j) = phi(P(i,j),Phi0(i,j));
            
        end
    end
end
end









