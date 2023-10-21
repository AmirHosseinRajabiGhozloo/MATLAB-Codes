function A = AavgFunction(BCN,BCE,BCS,BCW,dx,dy,dz,kx,ky,ndx,ndy)

Ax = zeros(ndx,ndy);
Ay = zeros(ndx,ndy);
A = zeros(ndx,ndy,4);

for j = 1:ndy
    for i = 1:ndx

        Ax(i,j) = dy(i,j)*dz(i,j);
        Ay(i,j) = dx(i,j)*dz(i,j);

    end
end

for j = 1:ndy
    for i = 1:ndx

        if BCN(i,j) == 100

            A(i,j,1) = 2*Ax(i,j)*Ax(i-1,j)*kx(i,j)*kx(i-1,j)/(Ax(i,j)*kx(i,j)*dx(i-1,j)+Ax(i-1,j)*kx(i-1,j)*dx(i,j));

        else

            A(i,j,1) = Ax(i,j)*kx(i,j)/dx(i,j);

        end

        if BCE(i,j) == 100

            A(i,j,2) = 2*Ay(i,j)*Ay(i,j+1)*ky(i,j)*ky(i,j+1)/(Ay(i,j)*ky(i,j)*dy(i,j+1)+Ay(i,j+1)*ky(i,j+1)*dy(i,j));

        else

            A(i,j,2) = Ay(i,j)*ky(i,j)/dy(i,j);

        end

        if BCS(i,j) == 100

            A(i,j,3) = 2*Ax(i,j)*Ax(i+1,j)*kx(i,j)*kx(i+1,j)/(Ax(i,j)*kx(i,j)*dx(i+1,j)+Ax(i+1,j)*kx(i+1,j)*dx(i,j));

        else

            A(i,j,3) = Ax(i,j)*kx(i,j)/dx(i,j);

        end

        if BCW(i,j) == 100

            A(i,j,4) = 2*Ay(i,j)*Ay(i,j-1)*ky(i,j)*ky(i,j-1)/(Ay(i,j)*ky(i,j)*dy(i,j-1)+Ay(i,j-1)*ky(i,j-1)*dy(i,j));

        else

            A(i,j,4) = Ay(i,j)*ky(i,j)/dy(i,j);

        end
    end
end
end