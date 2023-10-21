function zf=rooT_fcn(A,B,U,W)

Zt=[1 -(1+B-U) A-U*B-U-W^2 -A*B+B*(W^2)+W^2];
zf=roots(Zt);
kk=0;
for i=1:numel(zf)
    kk=kk+1;
    if isreal(zf(kk))==0
            zf(kk)=[]; 
            kk=kk-1;
    end
end
[row]=find(zf>B);
zf=zf(row);

% zfl=min(zf);
% zfv=max(zf);

% a1=zz(2);
% a2=zz(3);
% a3=zz(4);
% Q=(3*a2-a1^2)/9;
% J=(9*a1*a2-27*a3-2*a1^3)/54;
% D=Q^3+J^2;
% if D>0
%     M=J-sqrt(D);
%     N=J+sqrt(D);
%     if M<0 && N<0
%         z1=-((abs(N))^(1/3))-((abs(M))^(1/3))-a1/3;
%     elseif M>0 && N<0
%         z1=-((abs(N))^(1/3))+((abs(M))^(1/3))-a1/3;
%     elseif M<0 && N>0
%         z1=((abs(N))^(1/3))-((abs(M))^(1/3))-a1/3;
%     elseif M>0 && N>0
%         z1=((abs(N))^(1/3))+((abs(M))^(1/3))-a1/3;
%     end
%     z2=z1;
%     z3=z1;
% elseif D<0
%     tta=acosd(J/sqrt(-Q^3));
%     z1=2*sqrt(-Q)*cosd(tta/3)-a1/3;
%     z2=2*sqrt(-Q)*cosd(tta/3+120)-a1/3;
%     z3=2*sqrt(-Q)*cosd(tta/3+240)-a1/3;
% elseif D==0
%     z1=2*J^(1/3)-a1/3;
%     z2=-J^(1/3)-a1/3;
%     z3=z2;
% end
% z=[z1 z2 z3];

end

