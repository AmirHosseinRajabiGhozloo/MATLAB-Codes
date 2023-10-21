clc;
clear;

%% Note : if RBC or LBC = 100 : boundary is located between two blocks and two blocks are interconnected
%% Note : if RBC or LBC < 100 : boundary is constant flow rate condition or dp/dx = c and the value of RBC or LBC is equal to c
%% Note : if RBC or LBC > 100 : boundary is constant pressure and at boundary we have P = C and the value of RBC or LBC is equal to C

%% ==========================================================================================================================================================
%% ++++++++++++++++++++++++++ INPUTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ==========================================================================================================================================================

Property = xlsread('Properties');
dx = xlsread('dx');                                                                                                  %% Block dimension in x direction 
dy = xlsread('dy');                                                                                                  %% Block dimension in y direction 
dz = xlsread('dz');                                                                                                  %% Block dimension in z direction 
BCW = xlsread('BCW');                                                                                                %% West Boundary Condition
BCE = xlsread('BCE');                                                                                                %% East Boundary Condition
BCS = xlsread('BCS');                                                                                                %% South Boundary Condition
BCN = xlsread('BCN');                                                                                                %% East Boundary Condition
Internal_Boundary = Property(1,1);                                                                                   %% Type of internal boundary
BHP_o(:,:,1) = xlsread('BHP_o');                                                                                     %% oil Well Bottom Hole Pressure
BHP_w(:,:,1) = xlsread('BHP_w');                                                                                     %% water Well Bottom Hole Pressure
qosc(:,:,1) = xlsread('Qo');                                                                                         %% oil Well Flow Rate At Standard Condition
qwsc(:,:,1) = xlsread('Qw');                                                                                         %% water Well Flow Rate At Standard Condition
rw_o = xlsread('rw_o');                                                                                              %% oil Wellbore Radius
rw_w = xlsread('rw_w');                                                                                              %% water Wellbore Radius
So = xlsread('Skin_o');                                                                                              %% oil well Skin Factor
Sw = xlsread('Skin_w');                                                                                              %% water well Skin Factor
dt = Property(1,2);                                                                                                  %% Time step
ndt = Property(1,3);                                                                                                 %% Final time of simulation
Poi = Property(1,4);                                                                                                 %% Initial pressure
Co = Property(1,5);                                                                                                  %% Oil cmpressibility
Cw = Property(1,6);                                                                                                  %% water cmpressibility
Cf = Property(1,7);                                                                                                  %% rock cmpressibility
Bo0 = Property(1,8);                                                                                                 %% Reference oil FVF
Bw0 = Property(1,9);                                                                                                 %% Reference oil FVF
Rho_o = Property(1,10);                                                                                              %% oil density
Rho_w = Property(1,11);                                                                                              %% water density
Kx = xlsread('Kx');                                                                                                  %% X direction Permeability
Ky = xlsread('Ky');                                                                                                  %% Y direction permeability
Phi0 = xlsread('Phi0');                                                                                              %% initial porosity
P0 = Property(1,12);                                                                                                 %% refrence pressure
o = Property(1,13);                                                                                                  %% OIL WELL internal boundary condition
w = Property(1,14);                                                                                                  %% WATER WELL internal boundary condition
gama_o = 0.00021584*Rho_o*32.17;
gama_w = 0.00021584*Rho_w*32.17;
a = size(dx);
ndx = a(1);                                                                                                          %% Number of gridblocks in x direction
ndy = a(2);                                                                                                          %% Number of gridblocks in y direction
Aavg = AavgFunction(BCN,BCE,BCS,BCW,dx,dy,dz,Kx,Ky,ndx,ndy);
Aavg_N = Aavg(:,:,1);
Aavg_E = Aavg(:,:,2);
Aavg_S = Aavg(:,:,3);
Aavg_W = Aavg(:,:,4);
req = 0.14*(dx.^2+dy.^2).^(0.5);                                                                                     %% Equivalent drainage radius for each grid block 
dt = 2*dt;
Po = zeros(ndx,ndy,ndt+1);
Pw = zeros(ndx,ndy,ndt+1);
Po(1:ndx,1:ndy,1) = Poi;
Qosc = zeros(ndx,ndy,ndt+1);
Qwsc = zeros(ndx,ndy,ndt+1);
Qosc(:,:,1) = qosc;
Qwsc(:,:,1) = qwsc;
sWr = input('enter residual water saturation: '); 
sOr = input('enter residual oil saturation: ');
swi = input('enter initial water saturation: ');                                                                     %% should be higher than Swr
t(1) = 0;
sw = zeros(ndx,ndy,ndt+1);
so = zeros(ndx,ndy,ndt+1);
sw(1:ndx,1:ndy,1) = swi;
so(1:ndx,1:ndy,1) = 1-swi;
Phi = zeros(ndx,ndy,ndt);
Pc = zeros(ndx,ndy,ndt);
dPc = zeros(ndx,ndy,ndt);
dBo = zeros(ndx,ndy,ndt);
dBw = zeros(ndx,ndy,ndt);
ToE = zeros(ndx,ndy,ndt);
ToW = zeros(ndx,ndy,ndt);
TwE = zeros(ndx,ndy,ndt);
TwW = zeros(ndx,ndy,ndt);
ToS = zeros(ndx,ndy,ndt);
ToN = zeros(ndx,ndy,ndt);
TwS = zeros(ndx,ndy,ndt);
TwN = zeros(ndx,ndy,ndt);
Cpoo = zeros(ndx,ndy,ndt);
Cswo = zeros(ndx,ndy,ndt);
Cpow = zeros(ndx,ndy,ndt);
Csww = zeros(ndx,ndy,ndt);
alfa = zeros(ndx,ndy,ndt);
miu_o = zeros(ndx,ndy,ndt);
miu_w = zeros(ndx,ndy,ndt);
Bo = zeros(ndx,ndy,ndt);                                                                                                                     %% oil formation volume factor
Bw = zeros(ndx,ndy,ndt);                                                                                                                     %% water formation volume factor
pc = PcFunction(sWr,sOr,ndx,ndy,sw(:,:,1));
Pc(:,:,1) = pc(:,:,1);
Pw(:,:,1) = Po(:,:,1)-Pc(:,:,1);
miu = ViscosityFunction(Po(:,:,1),Pw(:,:,1),ndx,ndy);
miu_o(:,:,1) = miu(:,:,1);
miu_w(:,:,1) = miu(:,:,2);
B = FVFFunction(Po(:,:,1),Pw(:,:,1),ndx,ndy,Bo0,Bw0,Co,Cw,P0);                                                                                %% formation volume factor function
Bo(:,:,1) = B(:,:,1);                                                                                                                         %% oil formation volume factor
Bw(:,:,1) = B(:,:,2);                                                                                                                         %% water formation volume factor
Kr = KrFunction(sWr,sOr,ndx,ndy,sw(:,:,1));                                                                                                   %% relative permeability function
Kro(:,:,1) = Kr(:,:,1);                                                                                                                       %% oil relative permeability
Krw(:,:,1) = Kr(:,:,2);                                                                                                                       %% water relative permeability
landa = MobilityFunction(BCW,BCE,BCN,BCS,Kro(:,:,1),Krw(:,:,1),miu_w(:,:,1),miu_o(:,:,1),Bo(:,:,1),Bw(:,:,1),ndx,ndy,Po(:,:,1),Pw(:,:,1));    %% mobility function
landa_o(:,:,1) = landa(:,:,9);
landa_w(:,:,1) = landa(:,:,10);
Jo(:,:,1) = 2*pi*1.127*Kx.*dz.*landa_o(:,:,1)./(log(req./rw_o)+So);                                                                            %% oil wells Productivity index
Jw(:,:,1) = 2*pi*1.127.*Kx.*dz.*(Bo(:,:,1)./Bw(:,:,1).*landa_o(:,:,1)+landa_w(:,:,1))./(log(req./rw_w)+Sw);                                    %% water wells Productivity index
Jwo(:,:,1) = 2*pi*1.127.*Kx.*dz.*landa_w(:,:,1)./(log(req./rw_o)+So);                                                                          %% water productivity index at oil well site

if o == 1 && w == 1

    BHP_w(:,:,1) = Pw(:,:,1)+Qwsc(:,:,1)./Jw(:,:,1);
    BHP_o(:,:,1) = Po(:,:,1)+Qosc(:,:,1)./Jo(:,:,1);
    qwo(:,:,1) = Qosc(:,:,1).*landa_w(:,:,1)./landa_o(:,:,1);

for i = 1:ndx
    for j = 1:ndy

        if Qosc(i,j,1) == 0 

            BHP_o(i,j,1) = 0;

        end

        if Qwsc(i,j,1) == 0

            BHP_w(i,j,1) = 0;

        end
    end
end

elseif o == 1 && w == 2

    BHP_o(:,:,1) = Po(:,:,1)+Qosc(:,:,1)./Jo(:,:,1);
    Qwsc(:,:,1) = Jw(:,:,1).*(-Pw(:,:,1)+BHP_w(:,:,1));
    qwo(:,:,1) = Qosc(:,:,1).*landa_w(:,:,1)./landa_o(:,:,1);

for i = 1:ndx
    for j = 1:ndy

        if Qosc(i,j,1) == 0 

            BHP_o(i,j,1) = 0;

        end

        if BHP_w(i,j,1) == 0

            Qwsc(i,j,1) = 0;

        end
    end
end

elseif o == 2 && w == 1

    BHP_w(:,:,1) = Pw(:,:,1)+Qwsc(:,:,1)./Jw(:,:,1);
    Qosc(:,:,1) = -Jo(:,:,1).*(Po(:,:,1)-BHP_o(:,:,1));
    qwo(:,:,1) = -Jwo(:,:,1).*(Pw(:,:,1)-BHP_o(:,:,1));

for i = 1:ndx
    for j = 1:ndy

        if Qwsc(i,j,1) == 0
 
            BHP_w(i,j,1) = 0;

        end

        if BHP_o(i,j,1) == 0

            Qosc(i,j,1) = 0;

        end
    end
end

elseif o == 2 && w == 2

    Qosc(:,:,1) = -Jo(:,:,1).*(Po(:,:,1)-BHP_o(:,:,1));
    Qwsc(:,:,1) = Jw(:,:,1).*(-Pw(:,:,1)+BHP_w(:,:,1));
    qwo(:,:,1) = -Jwo(:,:,1).*(Pw(:,:,1)-BHP_o(:,:,1));

for i = 1:ndx
    for j = 1:ndy

        if BHP_w(i,j,1) == 0 

            Qwsc(i,j,1) = 0;

        end

        if BHP_o(i,j,1) == 0

            Qosc(i,j,1) = 0;

        end
    end
end
end

A = zeros(ndx*ndy,ndx*ndy,ndt);
b = zeros(ndx*ndy,ndt);

for n = 1:ndt

    t(n+1) = t(n)+dt;
    n
    Phi(:,:,n) = PorosityFunction(Po(:,:,n),Cf,Phi0,P0,ndx,ndy);
    pc = PcFunction(sWr,sOr,ndx,ndy,sw(:,:,n));
    miu = ViscosityFunction(Po(:,:,n),Pw(:,:,n),ndx,ndy);
    miu_o(:,:,n) = miu(:,:,1);
    miu_w(:,:,n) = miu(:,:,2);
    B = FVFFunction(Po(:,:,n),Pw(:,:,n),ndx,ndy,Bo0,Bw0,Co,Cw,P0);                                                                               %% formation volume factor function
    Bo(:,:,n) = B(:,:,1);                                                                                                                        %% oil formation volume factor
    Bw(:,:,n) = B(:,:,2);                                                                                                                        %% water formation volume factor
    Kr = KrFunction(sWr,sOr,ndx,ndy,sw(:,:,n));                                                                                                  %% relative permeability function
    Kro(:,:,n) = Kr(:,:,1);                                                                                                                      %% oil relative permeability
    Krw(:,:,n) = Kr(:,:,2);                                                                                                                      %% water relative permeability
    landa = MobilityFunction(BCW,BCE,BCN,BCS,Kro(:,:,n),Krw(:,:,n),miu_w(:,:,n),miu_o(:,:,n),Bo(:,:,n),Bw(:,:,n),ndx,ndy,Po(:,:,n),Pw(:,:,n));   %% mobility function
    dPc(:,:,n) = pc(:,:,2);
    dBo(:,:,n) = B(:,:,3);                                                                                                                       %% inverse oil formation volume factor derivative with respect to oil pressure
    dBw(:,:,n) = B(:,:,4);                                                                                                                       %% inverse water formation volume factor derivative with respect to water pressure
    landa_o_W = landa(:,:,1);                                                                                                                    %% oil mobility at west boundary for each grid
    landa_o_E = landa(:,:,2);                                                                                                                    %% oil mobility at east boundary for each grid
    landa_w_W = landa(:,:,5);                                                                                                                    %% water mobility at west boundary for each grid
    landa_w_E = landa(:,:,6);                                                                                                                    %% water mobility at east boundary for each grid
    landa_o_S = landa(:,:,4);                                                                                                                    %% oil mobility at west boundary for each grid    
    landa_o_N = landa(:,:,3);                                                                                                                    %% oil mobility at east boundary for each grid
    landa_w_S = landa(:,:,8);                                                                                                                    %% water mobility at west boundary for each grid
    landa_w_N = landa(:,:,7);                                                                                                                    %% water mobility at east boundary for each grid
    ToE(:,:,n) = 1.127*Aavg_E.*landa_o_E;                                                                                                        %% oil Transmissibility at East boundary
    ToW(:,:,n) = 1.127*Aavg_W.*landa_o_W;                                                                                                        %% oil Transmissibility at West boundary
    TwE(:,:,n) = 1.127*Aavg_E.*landa_w_E;                                                                                                        %% water Transmissibility at East boundary
    TwW(:,:,n) = 1.127*Aavg_W.*landa_w_W;                                                                                                        %% water Transmissibility at West boundary
    ToS(:,:,n) = 1.127*Aavg_S.*landa_o_S;                                                                                                        %% oil Transmissibility at South boundary
    ToN(:,:,n) = 1.127*Aavg_N.*landa_o_N;                                                                                                        %% oil Transmissibility at North boundary
    TwS(:,:,n) = 1.127*Aavg_S.*landa_w_S;                                                                                                        %% water Transmissibility at South boundary
    TwN(:,:,n) = 1.127*Aavg_N.*landa_w_N;                                                                                                        %% water Transmissibility at North boundary
    Cpoo(:,:,n) = dx.*dy.*dz.*(1-sw(:,:,n))/(5.615*dt).*(Cf*Phi0.*Bo(:,:,n).^(-1)+dBo(:,:,n).*Phi(:,:,n)); 
    Cswo(:,:,n) = -dx.*dy.*dz.*Phi(:,:,n)./(5.615*Bo(:,:,n)*dt);
    Cpow(:,:,n) = dx.*dy.*dz.*sw(:,:,n)/(5.615*dt).*(Cf*Phi0.*Bo(:,:,n).^(-1)+dBw(:,:,n).*Phi(:,:,n));
    Csww(:,:,n) = dx.*dy.*dz.*Phi(:,:,n)./(5.615*Bw(:,:,n)*dt)-dPc(:,:,n).*Cpow(:,:,n);
    alfa(:,:,n) = Cswo(:,:,n)./Csww(:,:,n);

    for i = 1:ndx
        for j=1:ndy

          if BCW(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCN(i,j) == 100

                if o == 1 && w == 1
                    
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n)));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n);
               
                elseif o == 1 && w == 2
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n)+Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
                
                elseif o == 2 && w == 1
                    
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
                
                elseif o == 2 && w == 2
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToS(i,j,n)+ToN(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n)+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n))+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i+1,j,n);
                
                end
                
                A(ndy*(i-1)+j,ndy*(i-1)+j-1,n) = ToW(i,j,n)-alfa(i,j,n)*TwW(i,j,n);
                A(ndy*(i-1)+j,ndy*(i-1)+j+1,n) = ToE(i,j,n)-alfa(i,j,n)*TwE(i,j,n);
                A(ndy*(i-1)+j,ndy*(i-1)+j+ndy,n) = ToS(i,j,n)-alfa(i,j,n)*TwS(i,j,n);
                A(ndy*(i-1)+j,ndy*(i-1)+j-ndy,n) = ToN(i,j,n)-alfa(i,j,n)*TwN(i,j,n);
           
 elseif BCW(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCN(i,j) == 100
               
                TwW(i,j,n) = 0;
                ToW(i,j,n) = 0;
               
                if o == 1 && w == 1
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n)));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n);
               
                elseif o == 1 && w == 2
                  
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n)+Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
                
                elseif o == 2 && w == 1
                    
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
                
                elseif o == 2 && w == 2
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToS(i,j,n)+ToN(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n)+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n))+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i+1,j,n);
               
                end
                
                A(ndy*(i-1)+j,ndy*(i-1)+j+1,n) = ToE(i,j,n)-alfa(i,j,n)*TwE(i,j,n);
                A(ndy*(i-1)+j,ndy*(i-1)+j+ndy,n) = ToS(i,j,n)-alfa(i,j,n)*TwS(i,j,n);
                A(ndy*(i-1)+j,ndy*(i-1)+j-ndy,n) = ToN(i,j,n)-alfa(i,j,n)*TwN(i,j,n);
            
 elseif BCW(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCN(i,j) == 100
                
                TwE(i,j,n) = 0;
                ToE(i,j,n) = 0;
                
                if o == 1 && w == 1
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n)));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n);
                
                elseif o == 1 && w == 2
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n)+Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
                
                elseif o == 2 && w == 1
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
               
                elseif o == 2 && w == 2
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToS(i,j,n)+ToN(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n)+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n))+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i+1,j,n);
              
                end
                
                A(ndy*(i-1)+j,ndy*(i-1)+j-1,n) = ToW(i,j,n)-alfa(i,j,n)*TwW(i,j,n);
                A(ndy*(i-1)+j,ndy*(i-1)+j+ndy,n) = ToS(i,j,n)-alfa(i,j,n)*TwS(i,j,n);
                A(ndy*(i-1)+j,ndy*(i-1)+j-ndy,n) = ToN(i,j,n)-alfa(i,j,n)*TwN(i,j,n);
            
 elseif BCW(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCN(i,j) == 100
                
                TwS(i,j,n) = 0;
                ToS(i,j,n) = 0;
                
                if o == 1 && w == 1
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n)));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n);
               
                elseif o == 1 && w == 2
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n)+Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n);
               
                elseif o == 2 && w == 1
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n);
               
                elseif o == 2 && w == 2
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToS(i,j,n)+ToN(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n)+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n))+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i+1,j,n);
               
                end
                
                A(ndy*(i-1)+j,ndy*(i-1)+j-1,n) = ToW(i,j,n)-alfa(i,j,n)*TwW(i,j,n);
                A(ndy*(i-1)+j,ndy*(i-1)+j+1,n) = ToE(i,j,n)-alfa(i,j,n)*TwE(i,j,n);
                A(ndy*(i-1)+j,ndy*(i-1)+j-ndy,n) = ToN(i,j,n)-alfa(i,j,n)*TwN(i,j,n);
           
 elseif BCW(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCN(i,j) < 100
               
                TwN(i,j,n) = 0;
                ToN(i,j,n) = 0;
                
                if o == 1 && w == 1
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n)));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
               
                elseif o == 1 && w == 2
                  
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n)+Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
               
                elseif o == 2 && w == 1
                  
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
                
                elseif o == 2 && w == 2
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToS(i,j,n)+ToN(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n)+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n))+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
               
                end
                
                A(ndy*(i-1)+j,ndy*(i-1)+j-1,n) = ToW(i,j,n)-alfa(i,j,n)*TwW(i,j,n);
                A(ndy*(i-1)+j,ndy*(i-1)+j+1,n) = ToE(i,j,n)-alfa(i,j,n)*TwE(i,j,n);
                A(ndy*(i-1)+j,ndy*(i-1)+j+ndy,n) = ToS(i,j,n)-alfa(i,j,n)*TwS(i,j,n);
          
 elseif BCW(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCN(i,j) == 100
               
                TwW(i,j,n) = 0;
                TwS(i,j,n) = 0;
                ToW(i,j,n) = 0;
                ToS(i,j,n) = 0;
               
                if o == 1 && w == 1
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n)));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n);
               
                elseif o == 1 && w == 2
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n)+Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n);
               
                elseif o == 2 && w == 1
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n);
               
                elseif o == 2 && w == 2
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToS(i,j,n)+ToN(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n)+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n))+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i+1,j,n);
               
                end
                
                A(ndy*(i-1)+j,ndy*(i-1)+j+1,n) = ToE(i,j,n)-alfa(i,j,n)*TwE(i,j,n);
                A(ndy*(i-1)+j,ndy*(i-1)+j-ndy,n) = ToN(i,j,n)-alfa(i,j,n)*TwN(i,j,n);
            
 elseif BCW(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCN(i,j) < 100
                
                TwW(i,j,n) = 0;
                TwN(i,j,n) = 0;
                ToW(i,j,n) = 0;
                ToN(i,j,n) = 0;
               
                if o == 1 && w == 1
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n)));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
               
                elseif o == 1 && w == 2
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n)+Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
               
                elseif o == 2 && w == 1
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
                
                elseif o == 2 && w == 2
                    
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToS(i,j,n)+ToN(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n)+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n))+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)-alfa(i,j,n)*TwE(i,j,n)*Pc(i,j+1,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
                
                end
                
                A(ndy*(i-1)+j,ndy*(i-1)+j+1,n) = ToE(i,j,n)-alfa(i,j,n)*TwE(i,j,n);
                A(ndy*(i-1)+j,ndy*(i-1)+j+ndy,n) = ToS(i,j,n)-alfa(i,j,n)*TwS(i,j,n);
           
 elseif BCW(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCN(i,j) == 100
               
                TwE(i,j,n) = 0;
                TwS(i,j,n) = 0;
                ToE(i,j,n) = 0;
                ToS(i,j,n) = 0;
               
                if o == 1 && w == 1
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n)));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n);
                
                elseif o == 1 && w == 2
                    
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n)+Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n);
               
                elseif o == 2 && w == 1
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i-1,j,n);
                
                elseif o == 2 && w == 2
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToS(i,j,n)+ToN(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n)+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n))+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwN(i,j,n)*Pc(i+1,j,n);
                
                end
                
                A(ndy*(i-1)+j,ndy*(i-1)+j-1,n) = ToW(i,j,n)-alfa(i,j,n)*TwW(i,j,n);
                A(ndy*(i-1)+j,ndy*(i-1)+j-ndy,n) = ToN(i,j,n)-alfa(i,j,n)*TwN(i,j,n);
           
 elseif BCW(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCN(i,j) < 100
                
                TwE(i,j,n) = 0;
                TwN(i,j,n) = 0;
                ToE(i,j,n) = 0;
                ToN(i,j,n) = 0;
               
                if o == 1 && w == 1
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n)));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
                
                elseif o == 1 && w == 2
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n));
                    b(ndy*(i-1)+j,n) = -Qosc(i,j,n)+alfa(i,j,n)*(Qosc(i,j,n)*landa_w(i,j,n)/landa_o(i,j,n)+Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
               
                elseif o == 2 && w == 1
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToN(i,j,n)+ToS(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Qwsc(i,j,n)+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
                
                elseif o == 2 && w == 2
                   
                    A(ndy*(i-1)+j,ndy*(i-1)+j,n) = -(ToW(i,j,n)+ToE(i,j,n)+ToS(i,j,n)+ToN(i,j,n)+Cpoo(i,j,n)-alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n)+Cpow(i,j,n))-alfa(i,j,n)*Jw(i,j,n)+Jo(i,j,n)-alfa(i,j,n)*Jwo(i,j,n));
                    b(ndy*(i-1)+j,n) = -Jo(i,j,n)*BHP_o(i,j,n)+alfa(i,j,n)*(Jw(i,j,n)*(BHP_w(i,j,n)+Pc(i,j,n))+Jwo(i,j,n)*(Pc(i,j,n)+BHP_o(i,j,n)))+(alfa(i,j,n)*Cpow(i,j,n)-Cpoo(i,j,n))*Po(i,j,n)+alfa(i,j,n)*(TwE(i,j,n)+TwW(i,j,n)+TwN(i,j,n)+TwS(i,j,n))*Pc(i,j,n)-alfa(i,j,n)*TwW(i,j,n)*Pc(i,j-1,n)-alfa(i,j,n)*TwS(i,j,n)*Pc(i+1,j,n);
               
                end
                
                A(ndy*(i-1)+j,ndy*(i-1)+j-1,n) = ToW(i,j,n)-alfa(i,j,n)*TwW(i,j,n);
                A(ndy*(i-1)+j,ndy*(i-1)+j+ndy,n) = ToS(i,j,n)-alfa(i,j,n)*TwS(i,j,n);
            
 end
        end
        end
        
    y = A(:,:,n)\b(:,n);
    Po(:,:,n+1) = reshape(y,ndx,ndy);
    Pw(:,:,n+1) = Po(:,:,n+1)-Pc(:,:,n);
    
    for i = 1:ndx
        for j = 1:ndy
            
            if BCW(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCN(i,j) == 100 
               
                sw(i,j,n+1) = sw(i,j,n)+1/Csww(i,j,n)*(TwW(i,j,n)*Pw(i,j-1,n+1)+TwE(i,j,n)*Pw(i,j+1,n+1)+TwS(i,j,n)*Pw(i+1,j,n)+TwN(i,j,n)*Pw(i-1,j,n+1)-(TwW(i,j,n)+TwE(i,j,n)+TwS(i,j,n)+TwN(i,j,n))*Pw(i,j,n+1)+Qwsc(i,j,n)+qwo(i,j,n)-Cpow(i,j,n)*(Po(i,j,n+1)-Po(i,j,n)));
            
            elseif BCW(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCN(i,j) == 100 
                
                sw(i,j,n+1) = sw(i,j,n)+1/Csww(i,j,n)*(TwE(i,j,n)*Pw(i,j+1,n+1)+TwS(i,j,n)*Pw(i+1,j,n)+TwN(i,j,n)*Pw(i-1,j,n+1)-(TwE(i,j,n)+TwS(i,j,n)+TwN(i,j,n))*Pw(i,j,n+1)+Qwsc(i,j,n)+qwo(i,j,n)-Cpow(i,j,n)*(Po(i,j,n+1)-Po(i,j,n)));
            
            elseif BCW(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCN(i,j) == 100 
                
                sw(i,j,n+1) = sw(i,j,n)+1/Csww(i,j,n)*(TwW(i,j,n)*Pw(i,j-1,n+1)+TwS(i,j,n)*Pw(i+1,j,n)+TwN(i,j,n)*Pw(i-1,j,n+1)-(TwW(i,j,n)+TwS(i,j,n)+TwN(i,j,n))*Pw(i,j,n+1)+Qwsc(i,j,n)+qwo(i,j,n)-Cpow(i,j,n)*(Po(i,j,n+1)-Po(i,j,n)));
            
            elseif BCW(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCN(i,j) == 100
                
                sw(i,j,n+1) = sw(i,j,n)+1/Csww(i,j,n)*(TwW(i,j,n)*Pw(i,j-1,n+1)+TwE(i,j,n)*Pw(i,j+1,n+1)+TwN(i,j,n)*Pw(i-1,j,n+1)-(TwW(i,j,n)+TwE(i,j,n)+TwN(i,j,n))*Pw(i,j,n+1)+Qwsc(i,j,n)+qwo(i,j,n)-Cpow(i,j,n)*(Po(i,j,n+1)-Po(i,j,n)));
            
            elseif BCW(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCN(i,j) < 100
                
                sw(i,j,n+1) = sw(i,j,n)+1/Csww(i,j,n)*(TwW(i,j,n)*Pw(i,j-1,n+1)+TwE(i,j,n)*Pw(i,j+1,n+1)+TwS(i,j,n)*Pw(i+1,j,n)-(TwW(i,j,n)+TwE(i,j,n)+TwS(i,j,n))*Pw(i,j,n+1)+Qwsc(i,j,n)+qwo(i,j,n)-Cpow(i,j,n)*(Po(i,j,n+1)-Po(i,j,n)));
            
            elseif BCW(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCN(i,j) < 100
                
                sw(i,j,n+1) = sw(i,j,n)+1/Csww(i,j,n)*(TwE(i,j,n)*Pw(i,j+1,n+1)+TwS(i,j,n)*Pw(i+1,j,n)-(TwE(i,j,n)+TwS(i,j,n))*Pw(i,j,n+1)+Qwsc(i,j,n)+qwo(i,j,n)-Cpow(i,j,n)*(Po(i,j,n+1)-Po(i,j,n)));
           
            elseif BCW(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCN(i,j) == 100
               
                sw(i,j,n+1) = sw(i,j,n)+1/Csww(i,j,n)*(TwE(i,j,n)*Pw(i,j+1,n+1)+TwN(i,j,n)*Pw(i-1,j,n+1)-(TwE(i,j,n)+TwN(i,j,n))*Pw(i,j,n+1)+Qwsc(i,j,n)+qwo(i,j,n)-Cpow(i,j,n)*(Po(i,j,n+1)-Po(i,j,n)));
            
            elseif BCW(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCN(i,j) < 100
               
                sw(i,j,n+1) = sw(i,j,n)+1/Csww(i,j,n)*(TwW(i,j,n)*Pw(i,j-1,n+1)+TwS(i,j,n)*Pw(i+1,j,n)-(TwW(i,j,n)+TwS(i,j,n))*Pw(i,j,n+1)+Qwsc(i,j,n)+qwo(i,j,n)-Cpow(i,j,n)*(Po(i,j,n+1)-Po(i,j,n)));
           
            elseif BCW(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCN(i,j) == 100
               
                sw(i,j,n+1) = sw(i,j,n)+1/Csww(i,j,n)*(TwW(i,j,n)*Pw(i,j-1,n+1)+TwN(i,j,n)*Pw(i-1,j,n+1)-(TwW(i,j,n)+TwN(i,j,n))*Pw(i,j,n+1)+Qwsc(i,j,n)+qwo(i,j,n)-Cpow(i,j,n)*(Po(i,j,n+1)-Po(i,j,n)));
            
            end
        end
    end
    
    so(:,:,n+1) = 1-sw(:,:,n+1);
    Phi(:,:,n+1) = PorosityFunction(Po(:,:,n+1),Cf,Phi0,P0,ndx,ndy);
    pc = PcFunction(sWr,sOr,ndx,ndy,sw(:,:,n+1));
    Pc(:,:,n+1) = pc(:,:,1);
    Pw(:,:,n+1) = Po(:,:,n+1)-Pc(:,:,n+1);
    miu = ViscosityFunction(Po(:,:,n+1),Pw(:,:,n+1),ndx,ndy);
    miu_o(:,:,n+1) = miu(:,:,1);
    miu_w(:,:,n+1) = miu(:,:,2);
    B = FVFFunction(Po(:,:,n+1),Pw(:,:,n+1),ndx,ndy,Bo0,Bw0,Co,Cw,P0); %% formation volume factor function
    Bo(:,:,n+1) = B(:,:,1); %% oil formation volume factor
    Bw(:,:,n+1) = B(:,:,2); %% water formation volume factor
    Kr = KrFunction(sWr,sOr,ndx,ndy,sw(:,:,n+1)); %% relative permeability function
    Kro(:,:,n+1) = Kr(:,:,1); %% oil relative permeability
    Krw(:,:,n+1) = Kr(:,:,2); %% water relative permeability
    landa = MobilityFunction(BCW,BCE,BCN,BCS,Kro(:,:,n+1),Krw(:,:,n+1),miu_w(:,:,n+1),miu_o(:,:,n+1),Bo(:,:,n+1),Bw(:,:,n+1),ndx,ndy,Po(:,:,n+1),Pw(:,:,n+1)); %% mobility function
    landa_o(:,:,n+1) = landa(:,:,9);
    landa_w(:,:,n+1) = landa(:,:,10);
    Jo(:,:,n+1) = 2*pi*1.127.*Kx.*dz.*landa_o(:,:,n+1)./(log(req./rw_o)+So); %% oil wells Productivity index
    Jw(:,:,n+1) = 2*pi*1.127.*Kx.*dz.*(Bo(:,:,n+1)./Bw(:,:,n+1).*landa_o(:,:,n+1)+landa_w(:,:,n+1))./(log(req./rw_w)+Sw); %% water wells Productivity index
    Jwo(:,:,n+1) = 2*pi*1.127.*Kx.*dz.*landa_w(:,:,n+1)./(log(req./rw_o)+So); %% water productivity index at oil well site
   
    if o == 1 && w == 1
        
        Qosc(:,:,n+1) = Qosc(:,:,n);
        Qwsc(:,:,n+1) = Qwsc(:,:,n);
        BHP_w(:,:,n+1) = Pw(:,:,n+1)+Qwsc(:,:,n+1)./Jw(:,:,n+1);
        BHP_o(:,:,n+1) = Po(:,:,n+1)+Qosc(:,:,n+1)./Jo(:,:,n+1);
        qwo(:,:,n+1) = Qosc(:,:,n+1).*landa_w(:,:,n+1)./landa_o(:,:,n+1);
        
        for i = 1:ndx
            for j = 1:ndy
                
                if Qosc(i,j,n+1) == 0 
                   
                    BHP_o(i,j,n+1) = 0;
               
                end
                
                if Qwsc(i,j,n+1) == 0
                  
                    BHP_w(i,j,n+1) = 0;
                
                end
            end
        end
        
    elseif o == 1 && w == 2
       
        Qosc(:,:,n+1) = Qosc(:,:,n);
        BHP_w(:,:,n+1) = BHP_w(:,:,n);
        BHP_o(:,:,n+1) = Po(:,:,n+1)+Qosc(:,:,n+1)./Jo(:,:,n+1);
        Qwsc(:,:,n+1) = Jw(:,:,n+1).*(-Pw(:,:,n+1)+BHP_w(:,:,n+1));
        qwo(:,:,n+1) = Qosc(:,:,n+1)*landa_w(:,:,n+1)./landa_o(:,:,n+1);
       
        for i = 1:ndx
            for j = 1:ndy
               
                if Qosc(i,j,n+1) == 0 
                   
                    BHP_o(i,j,n+1) = 0;
                
                end
                
                if BHP_w(i,j,n+1) == 0
                   
                    Qwsc(i,j,n+1) = 0;
                
                end
            end
        end
        
    elseif o == 2 && w == 1
       
        Qwsc(:,:,n+1) = Qwsc(:,:,n);
        BHP_o(:,:,n+1) = BHP_o(:,:,n);
        BHP_w(:,:,n+1) = Pw(:,:,n+1)+Qwsc(:,:,n+1)./Jw(:,:,n+1);
        Qosc(:,:,n+1) = -Jo(:,:,n+1).*(Po(:,:,n+1)-BHP_o(:,:,n+1));
        qwo(:,:,n+1) = -Jwo(:,:,n+1).*(Po(:,:,n+1)-BHP_o(:,:,n+1));
        
        for i = 1:ndx
             for j = 1:ndy
                 
                 if Qwsc(i,j,n+1) == 0 
                    
                     BHP_w(i,j,n+1) = 0;
                
                 end
                 
                 if BHP_o(i,j,n+1) == 0
                    
                     Qosc(i,j,n+1) = 0;
                
                 end
             end
        end
        
    elseif o == 2 && w == 2
        
        BHP_o(:,:,n+1) = BHP_o(:,:,n);
        BHP_w(:,:,n+1) = BHP_w(:,:,n);
        Qosc(:,:,n+1) = -Jo(:,:,n+1).*(Po(:,:,n+1)-BHP_o(:,:,n+1));
        Qwsc(:,:,n+1) = Jw(:,:,n+1).*(-Pw(:,:,n+1)+BHP_w(:,:,n+1));
        qwo(:,:,n+1) = -Jwo(:,:,n+1).*(Po(:,:,n+1)-BHP_o(:,:,n+1));
       
        for i = 1:ndx
            for j = 1:ndy
                
                if BHP_w(i,j,n+1) == 0 
                    
                    Qwsc(i,j,n+1) = 0;
                
                end
                
                if BHP_o(i,j,n+1) == 0
                   
                    Qosc(i,j,n+1) = 0;
               
                end
            end
        end
    end
    
    for i = 1:ndx
        for j = 1:ndy
            
            if abs(sw(i,j,n+1)-sWr) < 10^(-5)
               
                qwo(i,j,n+1) = 2;
            
            end
        end
    end
    end
    
    surf(sw(:,:,n));