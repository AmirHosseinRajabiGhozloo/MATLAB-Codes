clc;
clear;

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++ HAGEDORN & BROWN ++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Rho_l = input ('please enter the density of liquid in (lb/ft^3) ====>  ');
Rho_g = input ('please enter the density of gas in (lb/ft^3) ====>  ');
Vsl = input ('please enter the velocity of liquid in (ft/s) ====>  ');
Vsg = input ('please enter the velocity of gas in (ft/s) ====>  ');
Miul = input ('please enter the viscosity of liquid in (cp) ====>  ');
Miug = input ('please enter the viscosity of gas in (cp) ====>  ');
Sigmal = input('please enter the surface tension of liquid in (dynes/cm) ====>  ');
P = input('please enter the absolute pressure in (psia) ===>  ');
din = input ('please enter the inner diameter of the pipe in (in) ====>  ');
Roughness = input('please enter the roughness of the pipe if it is given (else enter 0) ====>  ');
Vm = Vsg+Vsl;
d = din/12;
Miu_l = Miul/1488.16;                                                      % lb/ft.s
Miu_g = Miug/1488.16;                                                      % lb/ft.s
Sigma_l = Sigmal/453.632;                                                  % lbm/s^3
g = 32.174;
gc = 32.174;
Nlv = Vsl*(Rho_l/(g*Sigma_l))^0.25;
Nl = Miu_l*(1/(Rho_l*(Sigma_l)^3))^0.25;
Nd = d*((g*Rho_l)/Sigma_l)^0.5;
Ngv = Vsg*(Rho_l/(g*Sigma_l))^0.25;

B = (Ngv*(Nl)^0.38)/(Nd)^2.14;
Saay = (1.0886-69.9473*B+2334.349*(B)^2-12896.638*(B)^3)/(1-53.4401*B+1517.936*(B)^2-8419.8115*(B)^3);

CNl = (0.0019+0.032*Nl-0.6642*(Nl)^2+4.995*(Nl)^3)/(1-10.01*Nl+33.869*(Nl)^2+277.28*(Nl)^3);
H = (Nlv/(Ngv)^0.575)*((P/14.7)^0.1)*(CNl/Nd);
fl = (((0.0047+1123.32*H+429489.64*(H)^2)/(1+1097.1566*H+722153.95*(H)^2))^0.5)*Saay;

Rho_m = (Rho_l*fl)+(Rho_g*(1-fl));
Miu_m = (Miu_l^fl)*(Miu_g^(1-fl));
NRe_m = Reynolds(Rho_m,Vm,d,Miu_m);
Fm = Friction_Factor(NRe_m,Roughness);

dp1 = Rho_m+((Fm*Rho_m*(Vm^2))/(2*gc*d));                                  % psf/ft
dp = dp1/148.73;                                                           % psi/ft

k = ['+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'];
disp(k)
y = [' ===>> dp/dz = ',num2str(dp),' psi/ft. <<=== '];
disp(y);
t = [' ===>> liquid holdup (fl) = 1 - gas hodup (fg) = ',num2str(fl),' <<==='];
disp(t);
disp(k);




