function A = AavgFunction(dx,dy,dz,Nx,Ny,Kx,Ky,BCN,BCE,BCS,BCW);

for j = 1:Ny
   
    for i = 1:Nx
        
        Ax(i,j) = dy(i,j)*dz(i,j);
        Ay(i,j) = dx(i,j)*dz(i,j);
        
    end
end

for j = 1:Ny
    
    for i = 1:Nx
     
     if BCN(i,j) == 100
         
        % if Kx(i,j) == 0 && Kx(i-1,j) == 0
             
          %   A(i,j,1) = 0;
             
         %else
           
         A(i,j,1) = 2*Ax(i,j)*Ax(i-1,j)*Kx(i,j)*Kx(i-1,j)/(Ax(i,j)*Kx(i,j)*dx(i-1,j)+Ax(i-1,j)*Kx(i-1,j)*dx(i,j));
        
         %end
         
     else
         
         A(i,j,1) = Ax(i,j)*Kx(i,j)/dx(i,j);
        
     end
     
     if BCE(i,j) == 100
         
        % if Ky(i,j) == 0 && Ky(i,j+1) == 0
             
            % A(i,j,2) = 0;
             
         %else
     
         A(i,j,2) = 2*Ay(i,j)*Ay(i,j+1)*Ky(i,j)*Ky(i,j+1)/(Ay(i,j)*Ky(i,j)*dy(i,j+1)+Ay(i,j+1)*Ky(i,j+1)*dy(i,j));
        
        % end
         
      else
         
         A(i,j,2) = Ay(i,j)*Ky(i,j)/dy(i,j);
        
     end
     
     if BCS(i,j) == 100
         
         %if Kx(i,j) == 0 && Kx(i+1,j) == 0
             
            % A(i,j,3) = 0;
             
        % else
    
         A(i,j,3) = 2*Ax(i,j)*Ax(i+1,j)*Kx(i,j)*Kx(i+1,j)/(Ax(i,j)*Kx(i,j)*dx(i+1,j)+Ax(i+1,j)*Kx(i+1,j)*dx(i,j));
         
         %end
        
     else
         
         A(i,j,3) = Ax(i,j)*Kx(i,j)/dx(i,j);
        
     end
     
     if BCW(i,j) == 100
         
        % if Ky(i,j) == 0 && Ky(i,j-1) == 0
             
            % A(i,j,4) = 0;
             
        % else
     
         A(i,j,4) = 2*Ay(i,j)*Ay(i,j-1)*Ky(i,j)*Ky(i,j-1)/(Ay(i,j)*Ky(i,j)*dy(i,j-1)+Ay(i,j-1)*Ky(i,j-1)*dy(i,j));
         
        % end
        
     else
         
         A(i,j,4) = Ay(i,j)*Ky(i,j)/dy(i,j);
        
     end           
    end
end
end

























