%% vicosity calculation
if state==1.5 || state==2.5 || state==0 % viscosity calculatin for gas phase
    Z=max(Zf);
    T=T*1.8;
    P=P/0.006895;
    Mw=sum(z.*(mw));
    X=3.5+(986/T)+0.01*Mw;
    Y=2.4-0.2*X;
    K=(9.4+0.02*Mw)*(T^1.5)/(2.9+19*Mw+T);
    dens=1.4935*0.001*(P*Mw)/(Z*T);
    viscosity=(10^-4)*K*exp(X*(dens^Y));
    
elseif state==1 || state==2  || state==0  % viscosity caculation for liquid from Lohrenz-Bray-clarc equation
    T=T*1.8; %R
    P=P/0.006895; %paia
    
    for i=1:n
        Tc(i)=Tc(i)*1.8;
        Pc(i)=Pc(i)/0.006895;
        mw(i)=mw(i);
        Tr(i)=T/Tc(i);
        eps=5.4402*(Tc(i)^(1/6))/((mw(i)^0.5)*(Pc(i)^(2/3)));
        
        if Tr(i)<=1.5
            vis(i)=34*(10^-5)*(Tr(i)^0.94)/eps;
        elseif Tr>1.5
            vis(i)=(17.78*(10^-5)*(4.58*Tr(i)-1.67)^0.625)/(eps);
        end
        
    end
    Mw=([component.Mw]);
    Ma=sum(zi.*Mw);
    Vc=([component.Vc])*(0.45359237/ 0.02831685)/Ma;
    ro_r=(sum(zi.*(Mw).*(Vc)))*(density)/Ma;
    vis_o=sum(zi.*(vis).*((Mw).^0.5))/sum(zi.*(Mw.^0.5));
    
    landa=(z'.*Vc)/(sum(z'.*Vc));
 Tcc=sum(landa.*Tc);
 TPc=sum(Tc.*z');
 PPc=sum(Pc.*z');
 wPc=sum(z'.*acf);
 Pcc=PPc*(1+(5.808+4.93*wPc)*((Tcc/TPc)-1));
 zPc=sum(z'.*Zc);
 Vcc=0.307*R*Tcc/Pcc
    Ppc=Ppc/0.006895;
    Tpc=Tpc*1.8;
    eps_m=5.4402*(Tpc^(1/6))/((Ma^0.5)*(Ppc^(2/3)));
    
    viscosity=vis_o+((0.1023+0.023364*ro_r+0.053533*(ro_r^2)-0.040758*(ro_r^3)+...
        0.0093324*(ro_r^4))^4-0.0001)/eps_m;
else
    disp('mixture is 3-phase and can not calculate property ')
    density=0;
    viscosity=0;
    
end