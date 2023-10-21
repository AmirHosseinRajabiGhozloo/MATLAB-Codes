function landa = MobilityFunction(BCW,BCE,BCN,BCS,Kro,Krw,miow,mioo,Bo,Bw,ndx,ndy,Po,Pw)

landa = zeros(ndx,ndy,10);

for i = 1:ndx
    for j = 1:ndy

        if BCW(i,j) == 100

            if Po(i,j) > Po(i,j-1)

                landa(i,j,1) = Kro(i,j)/(mioo(i,j)*Bo(i,j));

            else

                landa(i,j,1) = Kro(i,j-1)/(mioo(i,j-1)*Bo(i,j-1));

            end

            if Pw(i,j) > Pw(i,j-1)

                landa(i,j,5) = Krw(i,j)/(miow(i,j)*Bw(i,j));

            else

                landa(i,j,5) = Krw(i,j-1)/(miow(i,j-1)*Bw(i,j-1));

            end

        else

            landa(i,j,1) = Kro(i,j)/(mioo(i,j)*Bo(i,j));
            landa(i,j,5) = Krw(i,j)/(miow(i,j)*Bw(i,j));

        end

        if BCE(i,j) == 100

            if Po(i,j) > Po(i,j+1)

                landa(i,j,2) = Kro(i,j)/(mioo(i,j)*Bo(i,j));

            else

                landa(i,j,2) = Kro(i,j+1)/(mioo(i,j+1)*Bo(i,j+1));

            end

            if Pw(i,j) > Pw(i,j+1)

                landa(i,j,6) = Krw(i,j)/(miow(i,j)*Bw(i,j));

            else

                landa(i,j,6) = Krw(i,j+1)/(miow(i,j+1)*Bw(i,j+1));

            end

        else

            landa(i,j,2) = Kro(i,j)/(mioo(i,j)*Bo(i,j));
            landa(i,j,6) = Krw(i,j)/(miow(i,j)*Bw(i,j));

        end

        if BCN(i,j) == 100

            if Po(i,j) > Po(i-1,j)

                landa(i,j,3) = Kro(i,j)/(mioo(i,j)*Bo(i,j));

            else

                landa(i,j,3) = Kro(i-1,j)/(mioo(i-1,j)*Bo(i-1,j));

            end

            if Pw(i,j) > Pw(i-1,j)

                landa(i,j,7) = Krw(i,j)/(miow(i,j)*Bw(i,j));

            else

                landa(i,j,7) = Krw(i-1,j)/(miow(i-1,j)*Bw(i-1,j));

            end

        else

            landa(i,j,3) = Kro(i,j)/(mioo(i,j)*Bo(i,j));
            landa(i,j,7) = Krw(i,j)/(miow(i,j)*Bw(i,j));

        end

        if BCS(i,j) == 100

            if Po(i,j) > Po(i+1,j)

                landa(i,j,4) = Kro(i,j)/(mioo(i,j)*Bo(i,j));

            else

                landa(i,j,4) = Kro(i+1,j)/(mioo(i+1,j)*Bo(i+1,j));

            end

            if Pw(i,j) > Pw(i+1,j)

                landa(i,j,8) = Krw(i,j)/(miow(i,j)*Bw(i,j));

            else

                landa(i,j,8) = Krw(i+1,j)/(miow(i+1,j)*Bw(i+1,j));

            end

        else

            landa(i,j,4) = Kro(i,j)/(mioo(i,j)*Bo(i,j));
            landa(i,j,8) = Krw(i,j)/(miow(i,j)*Bw(i,j));

        end

        landa(i,j,9) = Kro(i,j)/(mioo(i,j)*Bo(i,j));
        landa(i,j,10) = Krw(i,j)/(miow(i,j)*Bw(i,j));

    end
end
end