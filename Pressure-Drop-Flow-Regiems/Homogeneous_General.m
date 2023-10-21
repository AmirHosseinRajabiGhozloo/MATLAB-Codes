clc;
clear;

% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++ HOMOGENEOUS FLOW +++++++++++++++++++++++++++++++++++++++
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
fg = Vsg/Vm;
fl = 1-fg;
Rho_m = 1/((x/Rho_g)+((1-x)/Rho_l));
Miu_m = Miu_l*(1-x)+Miu_g*(x);
NRe_m = Reynolds(Rho_m,Vm,d,Miu_m);
Fm = Friction_Factor(NRe_m,Roughness);

dp1 = Rho_m + ((Fm*Rho_m*Vm^2)/(2*gc*d));                                  % psf/ft
dp = dp1/148.73;                                                           % psi/ft
k = ['+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'];
disp(k)
y = [' ===>> dp/dz = ',num2str(dp),' psi/ft. <<=== '];
disp(y);
t = [' ===>> liquid holdup (fl) = 1 - gas hodup (fg) = ',num2str(fl),' <<==='];
disp(t);
disp(k);