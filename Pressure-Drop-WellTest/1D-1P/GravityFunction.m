function Qg = GravityFunction(N,Z,gama,LBC,RBC,LTx,RTx);

for i = 1:N
    
    if LBC(i) == 100
        
     LQg(i) = LTx(i)*gama*(Z(i-1)-Z(i)); 
     
    else
        
        LQg(i) = 0;
        
    end
    
    if  RBC(i) == 100
        
        RQg(i) = RTx(i)*gama*(Z(i+1)-Z(i)); 
        
    else 
        
        RQg(i) = 0;
        
    end
    
    Qg(i) = LQg(i)+RQg(i);
    
end
end
