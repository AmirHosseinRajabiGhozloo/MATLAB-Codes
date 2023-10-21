function Qg = GravityFunction(Gama,Z,Nx,Ny,TN,TE,TS,TW,BCN,BCE,BCS,BCW);

for j = 1:Ny    
    
    for i = 1:Nx
    
        if BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) == 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))-TN(i,j)*Gama*(Z(i,j)-Z(i-1,j))+TE(i,j)*Gama*(Z(i,j+1)-Z(i,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) > 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))-TN(i,j)*Gama*(Z(i,j)-Z(i-1,j))+TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) < 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))-TN(i,j)*Gama*(Z(i,j)-Z(i-1,j))+TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) == 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j))+TE(i,j)*Gama*(Z(i,j+1)-Z(i,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) == 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j))+TE(i,j)*Gama*(Z(i,j+1)-Z(i,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) == 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))-TN(i,j)*Gama*(Z(i,j)-Z(i-1,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) == 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))-TN(i,j)*Gama*(Z(i,j)-Z(i-1,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) == 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))+TE(i,j)*Gama*(Z(i,j+1)-Z(i,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) == 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))+TE(i,j)*Gama*(Z(i,j+1)-Z(i,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) == 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) == 100
        
            Qg(i,j) = TE(i,j)*Gama*(Z(i,j+1)-Z(i,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) > 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))+TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) == 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) > 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))-TN(i,j)*Gama*(Z(i,j)-Z(i-1,j));
        
        elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) > 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j))+TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) == 100
            
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) == 100
        
            Qg(i,j) = TE(i,j)*Gama*(Z(i,j+1)-Z(i,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) < 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))+TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) == 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) < 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))-TN(i,j)*Gama*(Z(i,j)-Z(i-1,j));
            
        elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) < 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j))+TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
        
        elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) == 100
            
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) == 100
        
            Qg(i,j) = TE(i,j)*Gama*(Z(i,j+1)-Z(i,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) > 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))+TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) == 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) > 100
            
                Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))-TN(i,j)*Gama*(Z(i,j)-Z(i-1,j));
        
        elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) > 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j))+TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) == 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) == 100
        
            Qg(i,j) = TE(i,j)*Gama*(Z(i,j+1)-Z(i,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) < 100
            
                Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))+TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
        
        elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) == 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j))-TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) < 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j))-TN(i,j)*Gama*(Z(i,j)-Z(i-1,j));
            
        elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) < 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j))+TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) > 100
        
            Qg(i,j)-TN(i,j)*Gama*(Z(i,j)-Z(i-1,j));
            
        elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) > 100
            
            Qg(i,j) = TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) > 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j));
            
        elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) == 100
        
            Qg(i,j) = -TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) < 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j));
            
        elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) < 100
        
            Qg(i,j) = TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) < 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j));
            
        elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) == 100
        
            Qg(i,j) = -TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) < 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j));
            
        elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) < 100
        
            Qg(i,j) = TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) < 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j));
            
        elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) == 100
        
            Qg(i,j) = -TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) > 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j));
            
        elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) > 100
        
            Qg(i,j) = TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) > 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j));
            
        elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) == 100
        
            Qg(i,j) = -TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) < 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j));
            
        elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) < 100
        
            Qg(i,j) = TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) < 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j));
            
        elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) == 100
        
            Qg(i,j) = -TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) > 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j));
            
        elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) > 100
        
            Qg(i,j) = TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) > 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j));
            
        elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) == 100
        
            Qg(i,j) = -TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) == 100
        
            Qg(i,j) = -TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) > 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j));
            
        elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) > 100
        
            Qg(i,j) = TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) > 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j));
            
        elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) == 100
        
            Qg(i,j) = -TW(i,j)*Gama*(Z(i,j)-Z(i,j-1));
            
        elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) < 100
        
            Qg(i,j) = TS(i,j)*Gama*(Z(i+1,j)-Z(i,j));
            
        elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) < 100
        
            Qg(i,j) = TE(i,j)*Gama*(Z(i,j+1)-Z(i,j));
            
        elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) < 100
        
            Qg(i,j) = -TN(i,j)*Gama*(Z(i,j)-Z(i-1,j));
            
        elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block
        
            Qg(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block
        
            Qg(i,j) = 0;
            
        elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block
        
            Qg(i,j) = 0;
            
        elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
        
            Qg(i,j) = 0;
            
        elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block
        
            Qg(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block
        
            Qg(i,j) = 0;
            
        elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block
        
            Qg(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block
        
            Qg(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block
        
            Qg(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
        
            Qg(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block
        
            Qg(i,j) = 0;
            
        elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block
        
            Qg(i,j) = 0;
            
        elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block
        
            Qg(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
        
            Qg(i,j) = 0;
            
        elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block
        
            Qg(i,j) = 0;
            
        elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
        
            Qg(i,j) = 0;
        
        end        
    end
end
end





