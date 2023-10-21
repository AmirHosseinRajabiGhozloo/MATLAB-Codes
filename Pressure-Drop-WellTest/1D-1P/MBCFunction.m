function MB = MBCFunction(dx,dy,dz,dt,Ngrid,LBC,RBC,L_Aavg,R_Aavg,P,P0,Co,B0,Phi0,Cphi,j,Qsc,Qqsc,initial_BC); 

if initial_BC == 1

    Qt = Qsc;
    
else
    
    Qt = Qqsc;                                                             
    
end

for k = 2:j+1
    
    Phi_a = PorosityFunction(P(:,k),Ngrid,P0,Phi0,Cphi);
    Phi_b = PorosityFunction(P(:,k-1),Ngrid,P0,Phi0,Cphi);
    Phi_c = PorosityFunction(P(:,1),Ngrid,P0,Phi0,Cphi);
    miu = ViscosityFunction(P(:,k),Ngrid,LBC,RBC);
    Lmiu = miu(:,1);                                                       %% viscosity at Left boundary
    Rmiu = miu(:,2);                                                       %% viscosity at Right boundary
    Bo_a = FVFFunction(Ngrid,P(:,k),Co,B0,P0,LBC,RBC);
    Bo_b = FVFFunction(Ngrid,P(:,k-1),Co,B0,P0,LBC,RBC);
    Bo_c = FVFFunction(Ngrid,P(:,1),Co,B0,P0,LBC,RBC);
    LBo = Bo_a(:,1);                                                       %% FVF at Left boundary
    RBo = Bo_a(:,2);                                                       %% FVF at Right boundary
    LT(:,k) = 1.127*L_Aavg./(Lmiu.*LBo);                                   %% Transmissibility at Left boundary
    RT(:,k) = 1.127*R_Aavg./(Rmiu.*RBo);                                   %% Transmissibility at Right boundary
    
    for i = 1:Ngrid
        
        if LBC(i) < 100
        
           LQ1(i,k) = (-LT(i,k)*LBC(i)*dx(i));                     
            
        elseif LBC(i) > 100
        
           LQ1(i,k) = (2*LT(i,k)*(LBC(i)-P(i,k)));
            
        else
            
           LQ1(i,k) = 0;
            
        end
    
        if RBC(i) < 100
        
           RQ1(i,k) = (RT(i,k)*RBC(i)*dx(i));                     
            
        elseif RBC(i) > 100
        
           RQ1(i,k) = (2*RT(i,k)*(RBC(i)-P(i,k)));
            
        else
            
           RQ1(i,k) = 0;
            
        end
        
           Q1(i,k) = ((LQ1(i,k)+RQ1(i,k)+Qt(i,k))*dt);
           Q2(i,k) = (dx(i)*dy(i)*dz(i)/5.615)*(Phi_a(i)/Bo_a(i,3)-Phi_b(i)/Bo_b(i,3)); 
    
    end
end

MB(1) = sum(Q2(:,j+1))/sum(Q1(:,j+1));
MB(2) = sum(sum(Q2))/sum(sum(Q1));
    
end   
 