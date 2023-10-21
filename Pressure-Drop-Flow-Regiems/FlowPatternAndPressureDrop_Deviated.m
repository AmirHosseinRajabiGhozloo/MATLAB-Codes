clc;
clear;

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++ DEVIATED WELLS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++ DISTINGUISHING FLOW PATTERN +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

aaa = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
bbb = ('===>>  It is a Deviated well or pipe.  <<===');
ccc = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
disp(aaa);
disp(bbb);
disp(ccc);

ddd = ('===>>  Which way do you prefer to distinguish flow pattern of this two phase flow ?  <<===');
disp(ddd);

n = input ('[[[ 1 - Beggs & Brill   or   2 - Mechanistic Method ]]] ===>  ');


disp(aaa);
eee = ('+++++++++++++ INPUTS +++++++++++++++++++++++++++++++++++++++');
disp(eee);
disp(aaa);

switch n
    
    case 1
        
        theta = input ('Please enter the Deviation of flow from horizontal in (degree) ===>  ')*pi/180;  % Radian
        P = input ('Please enter the pressure of the flow in (psig) ===>  ')+14.7;                       % psia
        T = input ('Please enter the temperature of the flow in (F) ===>  ')+460;                        % Rankin
        d = input ('Please enter the inner diameter of the pipe (in) ===>  ')/12;                        % ft
        Roughness = input('please enter the roughness of the pipe if it is given (else enter 0) ====>  ');
        Rho_o = input ('Please enter the density of oil in (lb/ft^3) ===>  ');                                           
        Rho_w = input ('Please enter the density of water in (lb/ft^3) ===>  ');
        Rho_g = input ('Please enter the density of gas in (lb/ft^3) ===>  ');
        qo = input ('Please enter the flow rate of oil in (STB/day) ===>  ');
        qw = input ('Please enter the flow rate of water in (STB/day) ===>  ');
        Sigma_o = input ('Please enter the surface tension of oil in (dyne/cm) ===>  ')*0.0022;           % lb/s^2
        Sigma_w = input ('Please enter the surface tension of water in (dyne/cm) ===>  ')*0.0022;         % lb/s^2
        Bo = input ('Please enter the oil formation volume factor in (bbl/STB) ===>  ');
        Bw = input ('Please enter the water formation volume factor in (bbl/STB) ===>  ');
        Miu_g = input ('Please enter the viscosity of gas in (cp) ===>  ')*0.000672;                      % lb/ft.s                           
        %Miu_w = input ('Please enter the viscosity of water in (cp) ===>  ')
        WOR = input ('Please enter the water oil ratio (WOR) ===>  ');
        GLR = input ('Please enter the gas liquid ratio (Rp) in (scf/STB) ===>  '); 
        Rs = input ('Please enter the solution gas ratio (Rs) in (scf/STBo) ===>  ');
        g = 32.14;                                                                                       % ft/s^2
        gc = 32.14;
        R = 10.732;
        Qo = qo*Bo*(5.615)/(24*3600);                                                                    % ft^3/s                                     
        Qw = qw*Bw*(5.615)/(24*3600);                                                                    % ft^3/s
        Ql = Qo+Qw;                                                                                      % ft^3/s
        fo = 1/(1+WOR);
        fw = 1-fo;
        Rho_l = Rho_o*fo+Rho_w*fw;
        Sigma_l = Sigma_o*fo+Sigma_w*fw;
        
        Miu_w = exp(1.003-(1.479*0.01*(T-460))+(1.982*0.00001*(T-460)^2))*0.000672;
        
        A1 = 0.3265; 
        A2 = -1.0700; 
        A3 = -0.5339; 
        A4 = 0.01569; 
        A5 = -0.05165; 
        A6 = 0.5475; 
        A7 = -0.7361; 
        A8 = 0.1844; 
        A9 = 0.1056; 
        A10 = 0.6134; 
        A11 = 0.7210;
        f111 = @(Gama_g) ((2.7*P*Gama_g)/(T*(Rho_g)))-(((0.27*(P/(756.8-131*Gama_g-3.6*(Gama_g^2)))/(((2.7*P*Gama_g)/(T*(Rho_g)))*(T/(169.2+349.5*Gama_g-74*(Gama_g^2)))))*(A1+(A2/(T/(169.2+349.5*Gama_g-74*(Gama_g^2))))+(A3/((T/(169.2+349.5*Gama_g-74*(Gama_g^2)))^3))+(A4/((T/(169.2+349.5*Gama_g-74*(Gama_g^2)))^4))+(A5/((T/(169.2+349.5*Gama_g-74*(Gama_g^2)))^5))))+(((0.27*(P/(756.8-131*Gama_g-3.6*(Gama_g^2)))/(((2.7*P*Gama_g)/(T*(Rho_g)))*(T/(169.2+349.5*Gama_g-74*(Gama_g^2)))))^2)*(A6+(A7/(T/(169.2+349.5*Gama_g-74*(Gama_g^2))))+(A8/((T/(169.2+349.5*Gama_g-74*(Gama_g^2)))^2))))-((A9*((0.27*(P/(756.8-131*Gama_g-3.6*(Gama_g^2)))/(((2.7*P*Gama_g)/(T*(Rho_g)))*(T/(169.2+349.5*Gama_g-74*(Gama_g^2)))))^5))*((A7/(T/(169.2+349.5*Gama_g-74*(Gama_g^2))))+(A8/((T/(169.2+349.5*Gama_g-74*(Gama_g^2)))^2))))+(A10*(1+A11*((0.27*(P/(756.8-131*Gama_g-3.6*(Gama_g^2)))/(((2.7*P*Gama_g)/(T*(Rho_g)))*(T/(169.2+349.5*Gama_g-74*(Gama_g^2)))))^2))*(((0.27*(P/(756.8-131*Gama_g-3.6*(Gama_g^2)))/(((2.7*P*Gama_g)/(T*(Rho_g)))*(T/(169.2+349.5*Gama_g-74*(Gama_g^2)))))^2)/((T/(169.2+349.5*Gama_g-74*(Gama_g^2)))^2))*(exp((-A11)*((0.27*(P/(756.8-131*Gama_g-3.6*(Gama_g^2)))/(((2.7*P*Gama_g)/(T*(Rho_g)))*(T/(169.2+349.5*Gama_g-74*(Gama_g^2)))))^2))))+1);                    
        Gama_g = NewtonRaphson(f111,0.5);
        Z = (2.7*P*Gama_g)/(T*Rho_g);
        Bg = 0.02827*Z*T/P;                                                                               % scf/STB
        Qg = Bg*qo*(GLR-Rs)/(24*3600);                                                                    % ft^3/s 
        
        Gama_o = Rho_o/Rho_w;
        API = (141.5/Gama_o)-131.5;  
        Miu_o = OilViscosity(T,P,API,Rs,Gama_g,Gama_o)*0.000672;                                          % lb/ft.s
        Miu_l = Miu_o*fo+Miu_w*fw;
        
        Ap = (pi/4)*(d^2);                                                                                % ft^2
        Vsl = Ql/Ap;                                                                                      % ft/s
        Vsg = Qg/Ap;                                                                                      % ft/s 
        Vm = Vsl+Vsg;                                                                                     % ft/s 
        
        Landa_l = Ql/(Ql+Qg);
        Landa_g = 1-Landa_l; 
        L1 = 316*(Landa_l^0.302);
        L2 = 0.0009252*(Landa_l^(-2.4684));
        L3 = 0.1*(Landa_l^(-1.4516));
        L4 = 0.5*(Landa_l^(-6.738));
        NFr = (Vm^2)/(g*d);
        
        if Landa_l < 0.01  &&  NFr < L1            
            
           aaaa = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           bbbb = ('===>>  It is a segregated flow pattern.  <<===');
           cccc = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           disp(aaaa);
           disp(bbbb);
           disp(cccc); 
           FP = 1;
           
        elseif   Landa_l >= 0.01  &&  NFr < L2
            
           aaaa = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           bbbb = ('===>>  It is a segregated flow pattern.  <<===');
           cccc = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           disp(aaaa);
           disp(bbbb);
           disp(cccc); 
           FP = 1;
           
        elseif   Landa_l < 0.4  &&  NFr >= L1
            
           aaaa = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           bbbb = ('===>>  It is a distributed flow pattern.  <<===');
           cccc = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           disp(aaaa);
           disp(bbbb);
           disp(cccc); 
           FP = 2;
           
        elseif   Landa_l >= 0.4  &&  NFr > L4
            
           aaaa = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           bbbb = ('===>>  It is a distributed flow pattern.  <<===');
           cccc = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           disp(aaaa);
           disp(bbbb);
           disp(cccc); 
           FP = 2;
           
        elseif   Landa_l >= 0.01 && Landa_l < 0.4  &&  NFr > L3 && NFr < L4
            
           aaaa = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           bbbb = ('===>>  It is an intermittent flow pattern.  <<===');
           cccc = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           disp(aaaa);
           disp(bbbb);
           disp(cccc);
           FP = 3;
           
        elseif   Landa_l >= 0.04  &&  NFr > L3 && NFr < L4
            
           aaaa = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           bbbb = ('===>>  It is an intermittent flow pattern.  <<===');
           cccc = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           disp(aaaa);
           disp(bbbb);
           disp(cccc);
           FP = 3;
           
        elseif   Landa_l >= 0.01  &&  NFr >= L2 && NFr <= L3
            
           aaaa = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           bbbb = ('===>>  It is a transition flow pattern.  <<===');
           cccc = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           disp(aaaa);
           disp(bbbb);
           disp(cccc);
           FP = 4;
           
        end
        
        if FP == 1
         
            a = 0.980; 
            b = 0.4846; 
            c = 0.0868; 
            dd = 0.011; 
            e = -3.768; 
            f = 3.539; 
            gg = -1.614;
            
        elseif FP == 2
            
            a = 1.065; 
            b = 0.5824; 
            c = 0.0609; 
            dd = 0; 
            e = 0; 
            f = 0; 
            gg = 0; 
            
        elseif FP == 3
           
            a = 0.845; 
            b = 0.5351; 
            c = 0.0173; 
            dd = 2.96; 
            e = 0.305; 
            f = -0.4473; 
            gg = 0.0978; 
            
        end
        
        HL_0 = (a*Landa_l^b)/(NFr^c); 
        NLv = 1.938*Vsl*(Rho_l/Sigma_l)^0.25; 
        C = (1-Landa_l)*log(dd*(Landa_l^e)*(NLv^f)*(NFr^gg)); 
        Saay = 1+C*(sin(1.8*theta)-1/3*(sin(1.8*theta))^3); 
        HL_theta = Saay*HL_0; 
        yy = Landa_l/(HL_theta^2);
        
        if 1 < yy  &&  yy < 1.2
            
            S = log(2.2*yy)-1.2;
            
        else
            
            S = log(yy)/(-0.0523+3.182*log(yy)-0.8725*((log(yy))^2)+0.01853*((log(yy))^4)); 
            
        end
        
        ftp_fn = exp(S);
        Rho_n = Rho_l*Landa_l+Rho_g*Landa_g;
        Miu_n = Miu_l*Landa_l+Miu_g*Landa_g;
        NRe_n = (Rho_n*Vm*d)/(Miu_n);
        fn = 1/(2*log10(NRe_n/(4.5223*log10(NRe_n)-3.8215)))^2;
        ftp = ftp_fn*fn;
        Rho_m = Rho_l*HL_theta+Rho_g*(1-HL_theta);
        dp = ((g/gc)*Rho_m*sin(theta))+ftp*Rho_n*(Vm^2)/(2*gc*d);                                        % Psi/m
        
        k = ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        disp(k)
        y = [' ===>> dp/dz = ',num2str(dp),' psi/ft. <<=== '];
        disp(y);
        t = [' ===>> liquid holdup (HL(theta)) =  ',num2str(HL_theta),' <<==='];
        disp(t);
        disp(k);

    
                
    case 2 
        
        theta = input ('Please enter the Deviation of flow from horizontal in (degree) ===>  ')*pi/180;  % Radian
        P = input ('Please enter the pressure of the flow in (psig) ===>  ')+14.7;                       % psia
        T = input ('Please enter the temperature of the flow in (F) ===>  ')+460;                        % Rankin
        d = input ('Please enter the inner diameter of the pipe (in) ===>  ')/12;                        % ft
        Roughness = input('please enter the roughness of the pipe if it is given (else enter 0) ====>  ');
        Rho_o = input ('Please enter the density of oil in (lb/ft^3) ===>  ');                                           
        Rho_w = input ('Please enter the density of water in (lb/ft^3) ===>  ');
        Rho_g = input ('Please enter the density of gas in (lb/ft^3) ===>  ');
        qo = input ('Please enter the flow rate of oil in (STB/day) ===>  ');
        qw = input ('Please enter the flow rate of water in (STB/day) ===>  ');
        Sigma_o = input ('Please enter the surface tension of oil in (dyne/cm) ===>  ')*0.0022;           % lb/s^2
        Sigma_w = input ('Please enter the surface tension of water in (dyne/cm) ===>  ')*0.0022;         % lb/s^2
        Bo = input ('Please enter the oil formation volume factor in (bbl/STB) ===>  ');
        Bw = input ('Please enter the water formation volume factor in (bbl/STB) ===>  ');
        Miu_g = input ('Please enter the viscosity of gas in (cp) ===>  ')*0.000672;                      % lb/ft.s                           
        %Miu_w = input ('Please enter the viscosity of water in (cp) ===>  ')
        WOR = input ('Please enter the water oil ratio (WOR) ===>  ');
        GLR = input ('Please enter the gas liquid ratio (Rp) in (scf/STB) ===>  '); 
        Rs = input ('Please enter the solution gas ratio (Rs) in (scf/STBo) ===>  ');
        g = 32.14;                                                                                       % ft/s^2
        gc = 32.14;
        R = 10.732;
        Qo = qo*Bo*(5.615)/(24*3600);                                                                    % ft^3/s                                     
        Qw = qw*Bw*(5.615)/(24*3600);                                                                    % ft^3/s
        Ql = Qo+Qw;                                                                                      % ft^3/s
        fo = 1/(1+WOR);
        fw = 1-fo;
        Rho_l = Rho_o*fo+Rho_w*fw;
        Sigma_l = Sigma_o*fo+Sigma_w*fw;
        
        Miu_w = exp(1.003-(1.479*0.01*(T-460))+(1.982*0.00001*(T-460)^2))*0.000672;                                                               % lb/ft.s
        
        A1 = 0.3265; 
        A2 = -1.0700; 
        A3 = -0.5339; 
        A4 = 0.01569; 
        A5 = -0.05165; 
        A6 = 0.5475; 
        A7 = -0.7361; 
        A8 = 0.1844; 
        A9 = 0.1056; 
        A10 = 0.6134; 
        A11 = 0.7210;
        f111 = @(Gama_g) ((2.7*P*Gama_g)/(T*(Rho_g)))-(((0.27*(P/(756.8-131*Gama_g-3.6*(Gama_g^2)))/(((2.7*P*Gama_g)/(T*(Rho_g)))*(T/(169.2+349.5*Gama_g-74*(Gama_g^2)))))*(A1+(A2/(T/(169.2+349.5*Gama_g-74*(Gama_g^2))))+(A3/((T/(169.2+349.5*Gama_g-74*(Gama_g^2)))^3))+(A4/((T/(169.2+349.5*Gama_g-74*(Gama_g^2)))^4))+(A5/((T/(169.2+349.5*Gama_g-74*(Gama_g^2)))^5))))+(((0.27*(P/(756.8-131*Gama_g-3.6*(Gama_g^2)))/(((2.7*P*Gama_g)/(T*(Rho_g)))*(T/(169.2+349.5*Gama_g-74*(Gama_g^2)))))^2)*(A6+(A7/(T/(169.2+349.5*Gama_g-74*(Gama_g^2))))+(A8/((T/(169.2+349.5*Gama_g-74*(Gama_g^2)))^2))))-((A9*((0.27*(P/(756.8-131*Gama_g-3.6*(Gama_g^2)))/(((2.7*P*Gama_g)/(T*(Rho_g)))*(T/(169.2+349.5*Gama_g-74*(Gama_g^2)))))^5))*((A7/(T/(169.2+349.5*Gama_g-74*(Gama_g^2))))+(A8/((T/(169.2+349.5*Gama_g-74*(Gama_g^2)))^2))))+(A10*(1+A11*((0.27*(P/(756.8-131*Gama_g-3.6*(Gama_g^2)))/(((2.7*P*Gama_g)/(T*(Rho_g)))*(T/(169.2+349.5*Gama_g-74*(Gama_g^2)))))^2))*(((0.27*(P/(756.8-131*Gama_g-3.6*(Gama_g^2)))/(((2.7*P*Gama_g)/(T*(Rho_g)))*(T/(169.2+349.5*Gama_g-74*(Gama_g^2)))))^2)/((T/(169.2+349.5*Gama_g-74*(Gama_g^2)))^2))*(exp((-A11)*((0.27*(P/(756.8-131*Gama_g-3.6*(Gama_g^2)))/(((2.7*P*Gama_g)/(T*(Rho_g)))*(T/(169.2+349.5*Gama_g-74*(Gama_g^2)))))^2))))+1);                    
        Gama_g = NewtonRaphson(f111,0.1);
        Z = (2.7*P*Gama_g)/(T*Rho_g);
        Bg = 0.02827*Z*T/P;                                                                               % scf/STB
        Qg = Bg*qo*(GLR-Rs)/(24*3600);                                                                    % ft^3/s 
        
        Gama_o = Rho_o/Rho_w;
        API = (141.5/Gama_o)-131.5;  
        Miu_o = OilViscosity(T,P,API,Rs,Gama_g,Gama_o)*0.000672;                                                 % lb/ft.s
        Miu_l = Miu_o*fo+Miu_w*fw;
        
        Ap = (pi/4)*(d^2);                                                                                % ft^2
        Vsl = Ql/Ap;                                                                                      % ft/s
        Vsg = Qg/Ap;                                                                                      % ft/s 
        Vm = Vsl+Vsg;                                                                                     % ft/s 
        
        Landa_l = Ql/(Ql+Qg);
        Landa_g = 1-Landa_l; 
        
        Vinf = 1.53*(g*Sigma_l*(Rho_l-Rho_g)/(Rho_l^2))^0.25;
        term1 = (0.43*Vsl+0.357*Vinf)*sin(theta);
        
        Vinft = (0.35*sin(theta)+0.54*cos(theta))*(g*d*(Rho_l-Rho_g)/Rho_l)^0.5;
        term2 = 12.19*(1.2*Vsl+Vinft);
        
        Eps = 0.2;
        
        if   Vsg < term1 
            
           aaaa = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           bbbb = ('===>>  It is a bubbly flow pattern.  <<===');
           cccc = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           disp(aaaa);
           disp(bbbb);
           disp(cccc);  
           
        elseif   Vsg >= term2-Eps  &&  Vsg <= term2+Eps
            
           aaaa = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           bbbb = ('===>>  It is a churn flow pattern.  <<===');
           cccc = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           disp(aaaa);
           disp(bbbb);
           disp(cccc); 
           
        else   
            
           aaaa = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           bbbb = ('===>>  It is a slug flow pattern.  <<===');
           cccc = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
           disp(aaaa);
           disp(bbbb);
           disp(cccc); 
           FP = 1;
            
        end 
        
        if FP == 1
            
            VinfT = (0.35*sin(theta)+0.54*cos(theta))*((g*d*(Rho_l-Rho_g))/(Rho_l))^0.5;
            
            if theta > (10*pi/180)  &&  theta <= (50*pi/180)
                
                C0 = 1.05;
                
            elseif theta > (50*pi/180)  &&  theta <= (60*pi/180)
                
                C0 = 1.15;
                
            elseif theta > (60*pi/180)  &&  theta <= (90*pi/180)
                
                C0 = 1.25;
                
            end
            
            VTB = C0*Vm+VinfT;
            fgLs = (Vsg)/(1.208*Vm+((1.41*((g*Sigma_l*(Rho_l-Rho_g))/(Rho_l^2))^0.25)*((sin(theta))^0.5)));
            flLs = 1-fgLs;
            VgLs = 1.08*Vm+((1.41*((g*Sigma_l*(Rho_l-Rho_g))/(Rho_l^2))^0.25)*((sin(theta))^0.5));
            VLLs = (Vm-(VgLs*fgLs))/(1-fgLs);
            ffff1 = @(flTB) (9.916*(g*d*(1-(1-flTB)^0.5))^0.5)-((((VTB-VLLs)*flLs)/(flTB))-VTB);
            flTB = NewtonRaphson(ffff1,0.1);
            fgTB = 1-flTB;
            VLTB = (9.916*(g*d*(1-(1-flTB)^0.5))^0.5);
            ffff11 = @(LLs_LSu) (VTB*LLs_LSu*(1-fgLs))+(VTB*(1-LLs_LSu)*(1-fgTB))-((VTB-VLLs)*(1-fgTB))-Vsl;
            LLs_LSu = NewtonRaphson(ffff11, 0.5);
            
            Miu_Ls = Miu_l*flLs+Miu_g*fgLs;
            Rho_Ls = Rho_l*flLs+Rho_g*fgLs;
            NRe_Ls = (Rho_Ls*Vm*d)/(Miu_Ls);
            fLs = Friction_Factor(NRe_Ls,Roughness);
            
            dp = (LLs_LSu*(Rho_l*flLs+Rho_g*fgLs)*((g*sin(theta))/(gc)))+(LLs_LSu*(fLs*Rho_Ls*(Vm^2))/(2*gc*d));
            
            k = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
            disp(k)
            y = [' ===>> dp/dz = ',num2str(dp),' psi/ft. <<=== '];
            disp(y);
            disp(k);
            
        end
        
end
