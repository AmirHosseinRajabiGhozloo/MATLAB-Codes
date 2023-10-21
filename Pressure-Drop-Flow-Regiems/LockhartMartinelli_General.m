clc;
clear;

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++ LOCKHART - MARTINELLI ++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Rho_l = input ('please enter the density of liquid in (lb/ft^3) ====>  ');
Rho_g = input ('please enter the density of gas in (lb/ft^3) ====>  ');
Vsl = input ('please enter the velocity of liquid in (ft/s) ====>  ');
Vsg = input ('please enter the velocity of gas in (ft/s) ====>  ');
Miul = input ('please enter the viscosity of liquid in (cp) ====>  ');
Miug = input ('please enter the viscosity of gas in (cp) ====>  ');
din = input ('please enter the inner diameter of the pipe in (in) ====>  ');
Roughness = input('please enter the roughness of the pipe if it is given (else enter 0) ====>  ');
Vm = Vsg+Vsl;
d = din/12;
Miu_l = Miul/1488.16;                                                      % lb/ft.s
Miu_g = Miug/1488.16;                                                      % lb/ft.s
x = (Vsg*Rho_g)/(Vsg*Rho_g+Vsl*Rho_l);
g = 32.174;
gc = 32.174;
X = ((((1-x)/x)^1.8)*(Rho_g/Rho_l)*((Miu_l/Miu_g)^0.2))^0.5;
fg = (1+(X)^0.8)^(-0.378);
fl = 1-fg;
NRe_g = Reynolds(Rho_g,Vsg,d,Miu_g);
NRe_l = Reynolds(Rho_l,Vsl,d,Miu_l);
Fl = Friction_Factor(NRe_l,Roughness);
Fg = Friction_Factor(NRe_g,Roughness);
dp_g = (Fg*Rho_g*(Vsg^2))/(2*gc*d);
dp_l = (Fl*Rho_l*(Vsl^2))/(2*gc*d);

if NRe_l >= 1000 & NRe_g >= 1000
    
    C = 20;
    
elseif NRe_l <= 1000 & NRe_g >= 1000
    
    C = 12;
    
elseif NRe_l >= 1000 & NRe_g <= 1000
    
    C = 10;
    
elseif NRe_l <= 1000 & NRe_g <= 1000
    
    C = 5;
     
end

Rho_m = Rho_g*fg+Rho_l*fl;
Phi_l = (1+(C/X)+(1/(X^2)))^0.5;
dp1 = Rho_m+(Phi_l^2)*dp_l;                                                % psf/ft
dp = dp1/148.73;                                                           % psi/ft

k = ['+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'];
disp(k)
y = [' ===>> dp/dz = ',num2str(dp),' psi/ft. <<=== '];
disp(y);
t = [' ===>> liquid holdup (fl) = 1 - gas hodup (fg) = ',num2str(fl),' <<==='];
disp(t);
disp(k);



