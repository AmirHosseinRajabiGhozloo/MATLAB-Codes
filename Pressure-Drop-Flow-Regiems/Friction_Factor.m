function f = Friction_Factor(NRe,Roughness)

if NRe <= 2100
    
    f = 64/NRe;
    
elseif NRe >= 2100
    
    if Roughness == 0
    
    f = 0.184*(NRe)^(-0.2);

    elseif Roughness > 0
    
    grad = (((Roughness)^(1.1098))/2.8257)+((7.149/NRe)^(0.8981));
    f = 4/((4*log10((Roughness/3.7065)-((5.0452/NRe)*log10(grad))))^(2));
    
    end
end

end