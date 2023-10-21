function MB = MBCFunction(P,dx,dy,dz,dt,Nx,Ny,n,IBC,Qsc,Qqsc,Co,B0,P0,Phi0,Cphi,N_Aavg,E_Aavg,S_Aavg,W_Aavg,BCN,BCE,BCS,BCW);

if IBC == 1
    
    Qt = Qsc;
    
else
    
    Qt = Qqsc;
    
end

for k = 2:n+1
    
   Phi_a = PorosityFunction(P(:,:,k),Nx,Ny,P0,Cphi,Phi0,BCN,BCE,BCS,BCW);
   Phi_b = PorosityFunction(P(:,:,k-1),Nx,Ny,P0,Cphi,Phi0,BCN,BCE,BCS,BCW);
   %Phi_c = PorosityFunction(P(:,:,1),Nx,Ny,P0,Cphi,Phi0,BCN,BCE,BCS,BCW);
   miu = ViscosityFunction(P(:,:,k),Nx,Ny,BCN,BCE,BCS,BCW);
   Nmiu = miu(:,:,1);
   Emiu = miu(:,:,2);
   Smiu = miu(:,:,3);
   Wmiu = miu(:,:,4);
   Bo_a = FVFFunction(P(:,:,k),Nx,Ny,Co,B0,P0,BCN,BCE,BCS,BCW);
   Bo_b = FVFFunction(P(:,:,k-1),Nx,Ny,Co,B0,P0,BCN,BCE,BCS,BCW);
   %Bo_c = FVFFunction(P(:,:,1),Nx,Ny,Co,B0,P0,BCN,BCE,BCS,BCW);
   NBo = Bo_a(:,:,1);
   EBo = Bo_a(:,:,2);
   SBo = Bo_a(:,:,3);
   WBo = Bo_a(:,:,4);
   TN(:,:,k) = 1.127*N_Aavg./(Nmiu.*NBo);
   TE(:,:,k) = 1.127*E_Aavg./(Emiu.*EBo);
   TS(:,:,k) = 1.127*S_Aavg./(Smiu.*SBo);
   TW(:,:,k) = 1.127*W_Aavg./(Wmiu.*WBo);
   
   for j = 1:Ny
       
       for i = 1:Nx
           
           if BCN(i,j) < 100
               
               Q1N(i,j,k) = -TN(i,j,k)*BCN(i,j)*dx(i,j);
            
           elseif BCN(i,j) > 100
           
               Q1N(i,j,k) = 2*TN(i,j,k)*(BCN(i,j)-P(i,j,k));
            
           else
               
               Q1N(i,j,k) = 0;
            
           end
           
           if BCE(i,j) < 100
           
               Q1E(i,j,k) = TE(i,j,k)*BCE(i,j)*dy(i,j);
            
           elseif BCE(i,j) > 100
           
               Q1E(i,j,k) = 2*TE(i,j,k)*(BCE(i,j)-P(i,j,k));
            
           else
               
               Q1E(i,j,k) = 0;
            
           end
           
           if BCS(i,j) < 100
           
               Q1S(i,j,k) = TS(i,j,k)*BCS(i,j)*dx(i,j);
            
           elseif BCS(i,j) > 100
           
               Q1S(i,j,k) = 2*TS(i,j,k)*(BCS(i,j)-P(i,j,k));
            
           else
               
               Q1S(i,j,k) = 0;
            
           end
           
           if BCW(i,j) < 100
           
               Q1W(i,j,k) = -TW(i,j,k)*BCW(i,j)*dy(i,j);                     
            
           elseif BCW(i,j) > 100
           
               Q1W(i,j,k) = 2*TW(i,j,k)*(BCW(i,j)-P(i,j,k));
            
           else
               
               Q1W(i,j,k) = 0;
            
           end
           
            Q1(i,j,k) = (Q1N(i,j,k)+Q1E(i,j,k)+Q1S(i,j,k)+Q1W(i,j,k)+Qt(i,j,k))*dt(n);
            Q2(i,j,k) = (dx(i,j)*dy(i,j)*dz(i,j)/5.615)*(Phi_a(i,j)/Bo_a(i,j,5)-Phi_b(i,j)/Bo_b(i,j,5));
            
            if BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block
               
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
            
            elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block
                
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
           
            elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block
                
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
            
            elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
               
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
           
            elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block
               
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
           
            elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block
              
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
           
            elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block
               
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
           
            elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block
               
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
           
            elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block
              
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
            
            elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
                
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
            
            elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block
               
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
           
            elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block
               
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
            
            elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block
                
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
            
            elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
               
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
           
            elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block
                
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
            
            elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
               
                Q1(i,j,k) = 0;
                Q2(i,j,k) = 0;
            
            end
       end
   end
end

t = ones(1,n+1)*dt(n);
MB(1) = sum(sum(Q2(:,:,n+1)))/sum(sum(Q1(:,:,n+1)));
MB(2) = sum(sum(sum(Q2)))/sum(sum(sum(Q1)));
            
end           
           
           
           
           
           
       
       
       
       
       
       
       
       
       
       
       
       








