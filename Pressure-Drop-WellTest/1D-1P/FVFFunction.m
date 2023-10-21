function  Bo = FVFFunction(N,P,C,B0,P0,LBC,RBC);

FVF = @(p)[B0/(1+C*(p-P0))];

for i = 1:N
    
    if LBC(i) == 100
      
        Bo(i,1) = FVF((P(i)+P(i-1))/2);
    
    else
        
        Bo(i,1) = FVF(P(i));
    
    end
    
     if RBC(i) == 100
        
        Bo(i,2) = FVF((P(i)+P(i+1))/2);
        
    else
        
        Bo(i,2) = FVF(P(i));
        
    end
    
    Bo(i,3) = FVF(P(i));

end
end
