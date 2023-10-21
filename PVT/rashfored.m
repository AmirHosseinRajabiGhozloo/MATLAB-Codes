function nv=Rashford_rice(zi,ki)


while true
    object=0;
syms sumr nv
sumr=abs(sum((zi.*(ki-1))./(1+(ki-1)*nv)));
left=1./(1-ki(1));
right=1./(1-ki(end));
start=(left+right)/2;
stop=0.5;
err=10^(-12);
while abs(start-stop)>err
    stop=start;
              start=vpa(start-(subs(sumr,nv,start))/(subs(diff(sumr),nv,start)));

    if double(start)>right 
       object=1;
        ki(1)=ki(1)*0.95;
         break
         
    else if  double(start)<left
            object=1;
         ki(1)=ki(1)*1.05;
         break    
    end
    
    end




end
if object==0
    break
end

end
clear nv
nv=double(vpa(start));
end