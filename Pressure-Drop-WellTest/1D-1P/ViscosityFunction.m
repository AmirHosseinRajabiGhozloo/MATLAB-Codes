function miu = ViscosityFunction(P,N,LBC,RBC);

viscosity=@(p)[10+2.5*10^(-5)*(p-14.7)];

for i = 1:N
    
    if LBC(i) == 100
        
 miu(i,1) = viscosity((P(i)+P(i-1))/2);
 
    else 
        
    miu(i,1) = viscosity(P(i));   
    
    end
    
    if RBC(i) == 100
        
    miu(i,2) = viscosity((P(i)+P(i+1))/2); 
    
    else
        
      miu(i,2) = viscosity(P(i));   
      
    end
    
    miu(i,3) = viscosity(P(i));  
    
end

end