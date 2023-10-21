function Miu_o = OilViscosity(T,P,API,Rs,Gama_g,Gama_o)

aa = 5.38088*10^(-3);
bb = 0.715082;
cc = -1.87784;
dd = 3.1437;
ee = 1.32657;
Pb = aa*(Rs^bb)*(Gama_g^cc)*(Gama_o^dd)*(T^ee);

A = 10^(3.0324-0.02023*API);
Miu_od = 10^(A*(T-460)^(-1.163))-1;

aa1 = 10.715*(Rs+100)^(-0.515);
bb1 = 5.44*(Rs+150)^(-0.338);
Miu_ob = (aa1*(Miu_od)^bb1);

AA = (-3.9*10^(-5))*P-5;
m = 2.6*(P^1.187)*(10^AA);
Miu_o = (Miu_ob*(P/Pb)^m);                                                 % cp

end