function A = AavgFunction(dx,dy,dz,Ngrid,Kx,LBC,RBC);

for i = 1:Ngrid
     
    Ax(i) = dy(i)*dz(i);

end

for i = 1:Ngrid
    
    if LBC(i) == 100
        
        A(i,1) = (2*Ax(i)*Ax(i-1)*Kx(i)*Kx(i-1))/(Ax(i)*Kx(i)*dx(i-1)+Ax(i-1)*Kx(i-1)*dx(i));
    
    else
        
        A(i,1) = (Ax(i)*Kx(i))/dx(i);
        
    end
    
    if RBC(i) == 100
        
        A(i,2) = (2*Ax(i)*Ax(i+1)*Kx(i)*Kx(i+1))/(Ax(i)*Kx(i)*dx(i+1)+Ax(i+1)*Kx(i+1)*dx(i));
        
    else
        
        A(i,2) = (Ax(i)*Kx(i))/dx(i);
        
    end
end
end

