or i=1:numel(w)
    Q(i)=0.25989-0.0217*acf(i)+0.00375*(acf(i)^2);
    zc(i)=1/(3*(1+Q(i)*acf(i)));
    if acf(i)<=0.4
        if ac(i)<=0.3671
            m0=0.465+1.347*acf(i)-0.528*acf(i)^2;
        else
            m0=0.5361+0.9593*acf(i);
        end
          
        m(i)=m0+0.01429*(5*(T/Tc(i))-3*m0-1)^2;
    elseif 0.4<acf(i)<0.55
        m0=0.5361+0.9593*acf(i);
        m1(i)=m0+0.01429*(5*(T/Tc(i))-3*m0-1)^2;
        m2(i)=m0+0.71*(T/Tc(i)-0.779)^2;
        m(i)=((0.55-acf(i))/0.15)*m1+((acf(i)-0.4)/0.5)*m2;
    else
        m(i)=m0+0.71*(T/Tc(i)-0.779)^2;
    end
    sayb(i)=zc(i)*Q(i);
    sayac(i)=(1-zc(i)*(1-Q(i)))^3;
    ac(i)=sayac(i)*(0.0083144*Tc(i)).^2./Pc(i);
    b(i)=sayb(i)*(0.0083144*Tc(i)/Pc(i));
    Tr(i)=T/Tc(i);
    if Tr(i)>1
        alpha(i)=1-(0.4774+1.328*w(i))*log(Tr(i));
    else
        alpha(i)=(1+m(i)*(1-Tr(i)^0.5))^2;
    end
    a(i)=alpha(i)*ac(i);
    
    
    
    am=0;
for i=1:n
    
    for j=1:n
        am=am+yi(i)*yi(j)*(((a(i)*a(j))^0.5));
    end
end

bm=sum(xi.*b);
B=bm*P/(0.0083144*T);
amm=sum(am);
A=amm*P/(0.0083144*T)^2;
%acfm=sum(x.*(b.^0.7).*w)/sum(x.*(b.^0.7));
u1=(1+3.*acf).*b;
u2=sum(u1.*xi);

w2=3.*acf.*(b.^2);
w1=w2.^0.5;
w1m=sum((w1).*xi);
W=(w1m)*P/(0.0083144*T);
U=u2*(P/(0.0083144*T));

zf=rooT_fcn(A,B,U,W);
zfv=max(zf);

for i = 1:n
    sumaij=0;
    for j = 1:n
        sumaij = sumaij + x(i) * x(j) * ((a(i) * a(j))^0.5);
    end
end
sumaij;

for i=1:n
    sumaij=0;
    for j=1:n
        sumaij=sumaij+(yi(j)*(((a(i)*a(j))^0.5)));
    end
     Phiv(i)=exp((b(i)/bm)*(zfv-1)-log(zfv-B)-(A/(B*(delTa2-delTa1)))*...
        (((2*sumaij)/am)-(b(i)/bm))*log((zfv+B*delTa2)/(zfv+B*delTa1)));
    ffv(i)=Phiv(i)*yi(i)*P;
end

 for i=1:n
    sumaij=0;
    for j=1:n
        sumaij = sumaij + x(j)*((a(i)*a(j))^0.5);
    end
p1(i)=-log(zfv-B) + (B* (b(i)/bm) ) / (zfv-B);
p2(i) =(A/(U^2 + 4*W^2)^0.5) .* ((2*sumaij /am)-(( (u1(i)/u2)*U^2 + 4*(w(i)/wm)*W^2 ) / (U^2 + 4*W^2))) .*log(  ( 2*zfv + U-sqrt(U^2+4*W^2) ) / (2*Z+U+sqrt(U^2+4*W^2) )  )  ;
p3(i)=A*(((2*(2*zfv+U)*(w(i)/w_mix)*W^2 )+ ((U*Z-2*W^2)*(uu(i)/u_mix )*U )) /(  (Z^2+U*Z-W^2)*(U^2+4*W^2)  )  )   ;
fug(i)=p1(i)+p2(i)-p3(i);
phi(i)=exp(fug(i)) ;
      
end
fug_coefficient=phi;
