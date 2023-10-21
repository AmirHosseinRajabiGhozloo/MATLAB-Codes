clc;
clear;

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++ VERTICAL WELLS ++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++ DISTINGUISHING FLOW PATTERN +++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

aaa = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
bbb = ('===>>  It is a Vertical well or pipe.  <<===');
ccc = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
disp(aaa);
disp(bbb);
disp(ccc);

Rho_l = input ('please enter the density of liquid in (Kg/m^3) ====>  ');
Rho_g = input ('please enter the density of gas in (Kg/m^3) ====>  ');
Vsl = input ('please enter the velocity of liquid in (m/s) ====>  ');
Vsg = input ('please enter the velocity of gas in (m/s) ====>  ');
Miu_l = input ('please enter the viscosity of liquid in (Pa.s) ====>  ');
Miu_g = input ('please enter the viscosity of gas in (Pa.s) ====>  ');
Sigma_l = input('please enter the surface tension of liquid in (N/m) ====>  ');
%P = input('please enter the absolute pressure in (psia) ===>  ');
din = input ('please enter the inner diameter of the pipe in (in) ====>  ');
Roughness = input('please enter the roughness of the pipe if it is given (else enter 0) ====>  ');
Vm = Vsg+Vsl;
d = din/39.3;
g = 9.8;                                                                   % m/s^2
gc = 9.8;

term1 = 3.1*((Sigma_l*g*(Rho_l-Rho_g))/(Rho_g^2))^0.25;

E = 1-exp((-0.125)*((10000*Vsg*Miu_g/Sigma_l)*((Rho_g/Rho_l)^0.5)-1.5));
Elc = (Vsl*E)/(Vsg+Vsl*E);
NRef = Reynolds(Rho_l,Vsl,d,Miu_l)*(1-E);
Ff = Friction_Factor(NRef,Roughness);
NResl = Reynolds(Rho_l,Vsl,d,Miu_l);
Fsl = Friction_Factor(NResl,Roughness);
dp_fl = (Fsl*Rho_l*(Vsl^2))/(2*d);
Vsc = Vsg+Vsl*E;
Rho_c = Rho_l*Elc+Rho_g*(1-Elc);
Miu_c = Miu_l*Elc+Miu_g*(1-Elc);
NRec = Reynolds(Rho_c,Vsc,d,Miu_c);
Fc = Friction_Factor(NRec,Roughness);
dp_c = (Fc*Rho_c*(Vsc^2))/(2*d);
Xm = (((1-E)^2)*(Ff/Fsl)*(dp_fl/dp_c))^0.5;
Ym = (g*(Rho_l-Rho_c))/(dp_c);

if E < 0.9
   
f1 = @(Delta_bar) (1+24*((Rho_l/Rho_g)^(1/3))*Delta_bar)/((((1-2*Delta_bar)^2)^2.5)*(1-(1-2*Delta_bar)^2))-(Ym+(Xm^2)/(1-(1-2*Delta_bar)^2)^3);
Delta_bar = NewtonRaphson(f1,0.2);
flf = 1-(1-2*Delta_bar)^2;

else
  
f1 = @(Delta_bar) (1+300*Delta_bar)/((((1-2*Delta_bar)^2)^2.5)*(1-(1-2*Delta_bar)^2))-(Ym+(Xm^2)/(1-(1-2*Delta_bar)^2)^3);
Delta_bar = NewtonRaphson(f1,0.2);
flf = 1-(1-2*Delta_bar)^2;    
    
end

term2 = flf+Elc*(1-flf);

fun = @(flf_min) (((1-E)^2)*(Ff/Fsl)*(Xm^2)*((2-1.5*flf_min)/((flf_min^3)*(1-1.5*flf_min)))) - Ym;
flf_min0 = 0.2;
options = optimset('Display','iter'); 
[flf_min,fval,exitflag,output] = fzero(fun,flf_min0,options);

NRevm = Reynolds(Rho_l,Vm,d,Miu_l);
Fvm = Friction_Factor(NRevm,Roughness);

f3 = @(Vm_calc) (2*(Vm_calc^1.2)*((Fvm/(2*d))^0.4)*((Rho_l/Sigma_l)^0.6)*(((0.4*Sigma_l)/(g*(Rho_l-Rho_g)))^0.5))-(0.725+4.15*(Vsg/Vm_calc)^0.5);
Vm_calc = NewtonRaphson(f3,0.2);

Vinf = 1.53*((g*(Rho_l-Rho_g)*Sigma_l)/(Rho_l^2))^0.25;
term3 = 0.429*Vsl+0.357*Vinf;

term4 = 19.1*((Sigma_l)/(g*(Rho_l-Rho_g)))^0.5;

if Vsg > term1  &&  term2 < 0.12  &&  flf < flfmin 
    
    aa = ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    bb = ('===>  It is an annular flow pattern.  <===');
    cc = ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    disp(aa);
    disp(bb);
    disp(cc);
    
elseif Vm_calc < Vm  &&  Vsg < 1.08*Vsl
    
    aa = ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    bb = ('===>  It is a dispersed bubbly flow pattern.  <===');
    cc = ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    disp(aa);
    disp(bb);
    disp(cc); 
    
elseif Vsg <= term3  &&  d > term4 
    
    aa = ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    bb = ('===>  It is a bubble flow pattern.  <===');
    cc = ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    disp(aa);
    disp(bb);
    disp(cc); 
    
else
    
    aa = ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    bb = ('===>  It is a slug flow pattern.  <===');
    cc = ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    disp(aa);
    disp(bb);
    disp(cc);
    FP = 1;
    
end

if FP == 1
   
    VinfT = 0.35*((g*d*(Rho_l-Rho_g))/(Rho_l))^0.5;
    VTB = 1.2*Vm+VinfT;
    fgLs = Vsg/(2.65*Vm+0.425);
    flLs = 1-fgLs;
    A_bar = VTB*fgLs+(1-fgLs)*(Vm-fgLs*(VinfT*flLs)^0.5);
    fff11 = @(fgTB) A_bar+((9.916*(g*d)^0.5)*(1-fgTB)*(1-(fgTB)^0.5))-(VTB*fgTB);
    fgTB = NewtonRaphson(fff11, 0.1);
    flTB = 1-fgTB;
    VlTB = 9.916*(g*d*(1-(fgTB)^0.5))^0.5;
    VLLs = VTB-((flTB*(VTB+VlTB))/(1-flTB));
    VgLs = 1.2*Vm+VinfT*(flLs^0.5);
    VgTB = VTB-(((1-flLs)*(1.2*Vm+VinfT-VgLs))/(1-flTB));
    ffff111 = @(Beta) Beta*VgTB*(1-flTB)+(1-Beta)*(1-flLs)*VgLs-Vsg;
    Beta = NewtonRaphson(ffff111, 0.1);
    
    Miu_Ls = Miu_l*flLs+Miu_g*fgLs;
    Rho_Ls = Rho_l*flLs+Rho_g*fgLs;
    NRe_Ls = Reynolds(Rho_Ls,Vm,d,Miu_Ls);
    fLs = Friction_Factor(NRe_Ls,Roughness);
    
    dp = ((g/gc)*(Rho_Ls*(1-Beta)+Beta*Rho_g))+((fLs*Rho_Ls*(1-Beta)*(Vm^2))/(2*gc*d));
    
    k = ('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
            disp(k)
            y = [' ===>> dp/dz = ',num2str(dp),' psi/m. <<=== '];
            disp(y);
            disp(k);
    
end










