clc;
clear;

%% Note : if RBC or LBC = 100 : boundary is located between two blocks and two blocks are interconnected
%% Note : if RBC or LBC < 100 : boundary is constant flow rate condition or dp/dx = c and the value of RBC or LBC is equal to c
%% Note : if RBC or LBC > 100 : boundary is constant pressure and at boundary we have P = C and the value of RBC or LBC is equal to C

%% ==========================================================================================================================================================
%% ++++++++++++++++++++++++++ INPUTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ==========================================================================================================================================================

Niteration = input('How many Iterations Do you Prefer to Run ? ');
Eps = input('How much Pressure diffrence is acceptable for breaking loop ? '); 
u = input('Please enter Pressure Update Factor ? ');
BHP(:,:,1) = xlsread('BottemHolePressure');
BCN = xlsread('NorthBoundryCondition');
BCE = xlsread('EastBoundryCondition');
BCS = xlsread('SouthBoundryCondition');
BCW = xlsread('WestBoundryCondition');
dx = xlsread('d_x');
dy = xlsread('d_y');
dz = xlsread('d_z');
Phi0 = xlsread('GridPorosity');
S = xlsread('Skin');
Top = xlsread('TopDepth');
Qsc = xlsread('WellQsc');
rw = xlsread('WellRadious');
Kx = xlsread('XPermeability');
Ky = xlsread('YPermeability');
inputs = xlsread('OtherProperties');
Co = inputs(1,1);
dt = inputs(1,2);
Tt = inputs(1,3);
Pi = inputs(1,4);
P0 = inputs(1,5);
IBC = inputs(1,6);
Rho = inputs(1,7);
Cphi = inputs(1,8);
B0 = inputs(1,9);
Z = Top + (dz/2);
b = size(dx);
Nx = b(1);
Ny = b(2);
Aavg = AavgFunction(dx,dy,dz,Nx,Ny,Kx,Ky,BCN,BCE,BCS,BCW);
N_Aavg = Aavg(:,:,1);
E_Aavg = Aavg(:,:,2);
S_Aavg = Aavg(:,:,3);
W_Aavg = Aavg(:,:,4);
Gama = (0.00021584*Rho*32.17);
Nt = (Tt/dt);
P(1:Nx,1:Ny,1) = Pi;
miu = ViscosityFunction(P(:,:,1),Nx,Ny,BCN,BCE,BCS,BCW);
Boi = FVFFunction(P(:,:,1),Nx,Ny,Co,B0,P0,BCN,BCE,BCS,BCW);


for i = 1:Ny
   
    for j = 1:Nx
       
        if BHP(i,j,1) == 0
       
            D(i,j) = 0;
       
        else
            
            D(i,j) = 1;
        
        end
    end
end

for j = 1:Ny
    
    for i = 1:Nx
        
if rw(i,j) > 0
   
    J(i,j) = (2*pi*1.127*(Kx(i,j)*Ky(i,j))^0.5*dz(i,j)/((log(0.2*dx(i,j)/rw(i,j))+S(i,j))*miu(i,j,5)*Boi(i,j,5))); 
  
else
    
    J(i,j) = 0;

end
    end
end

if IBC == 2
   
    Qqsc(:,:,1) = (-J.*(P(:,:,1)-BHP(:,:,1)));

end

t(1) = 0;

for n = 1:Nt
    
 Phi(:,:,n) = PorosityFunction(P(:,:,n),Nx,Ny,P0,Cphi,Phi0,BCN,BCE,BCS,BCW);
 miu = ViscosityFunction(P(:,:,n),Nx,Ny,BCN,BCE,BCS,BCW); 
 Nmiu(:,:,n) = miu(:,:,1);
 Emiu(:,:,n) = miu(:,:,2);
 Smiu(:,:,n) = miu(:,:,3);
 Wmiu(:,:,n) = miu(:,:,4);
 Bo = FVFFunction(P(:,:,n),Nx,Ny,Co,B0,P0,BCN,BCE,BCS,BCW);
 NBo(:,:,n) = Bo(:,:,1);
 EBo(:,:,n) = Bo(:,:,2);
 SBo(:,:,n) = Bo(:,:,3);
 WBo(:,:,n) = Bo(:,:,4);
 TN(:,:,n) = (1.127*N_Aavg./(Nmiu(:,:,n).*NBo(:,:,n)));
 TE(:,:,n) = (1.127*E_Aavg./(Emiu(:,:,n).*EBo(:,:,n)));
 TS(:,:,n) = (1.127*S_Aavg./(Smiu(:,:,n).*SBo(:,:,n)));
 TW(:,:,n) = (1.127*W_Aavg./(Wmiu(:,:,n).*WBo(:,:,n)));
 
 for j = 1:Ny
     
     for i = 1:Nx
         
         if rw(i,j) > 0
                      
             J(i,j) = (2*pi*1.127*(Kx(i,j)*Ky(i,j))^0.5*dz(i,j)/((log(0.2*dx(i,j)/rw(i,j))+S(i,j))*miu(i,j,5)*Bo(i,j,5)));
             
         else
             
             J(i,j) = 0;
         
         end
     end
 end
         
 Qg = GravityFunction(Gama,Z,Nx,Ny,TN(:,:,n),TE(:,:,n),TS(:,:,n),TW(:,:,n),BCN,BCE,BCS,BCW);
 Pg1 = P(:,:,n);
 
 for iteration = 1:Niteration
     
     Boi = FVFFunction(Pg1,Nx,Ny,Co,B0,P0,BCN,BCE,BCS,BCW);   
     Landa = (dx.*dy.*dz./(5.615).*(Phi(:,:,n)*Co/B0+Cphi*Phi0./Boi(:,:,5)));
     F = Landa./(TN(:,:,n)+TE(:,:,n)+TS(:,:,n)+TW(:,:,n));
     x = min(F);
     x = min(x);
     dt(n) = x;
     Cc = (dx.*dy.*dz./(5.615*dt(n)).*(Phi(:,:,n)*Co/B0+Cphi*Phi0./Boi(:,:,5))); 
     
     for j = 1:Ny
         
         for i = 1:Nx
             
             if BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));

             elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) > 100
 
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                         
             elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) < 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*P(i,j,n)*BCW(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                          
             elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) == 100
                
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
               
             elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                 
             elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*BCE(i,j)*dy(i,j)-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) == 100
              
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*BCN(i,j)*dx(i,j)+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) == 100
              
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));                            

             elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) > 100
               
               P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                          
             elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) > 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) > 100
                
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*BCN(i,j)*dx(i,j)+TE(i,j,n)*BCE(i,j)*dy(i,j)-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-TN(i,j,n)*BCN(i,j)*dx(i,j)+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                            
             elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) < 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*BCN(i,j)*dx(i,j)+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                            
             elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*BCE(i,j)*dy(i,j)-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) < 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*BCE(i,j)*dy(i,j)-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                          
             elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) < 100
               
               P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                          
             elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*BCN(i,j)*dx(i,j)+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                          
             elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-TN(i,j,n)*BCN(i,j)*dx(i,j)+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) > 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*BCN(i,j)*dx(i,j)+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*BCE(i,j)*dy(i,j)-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) > 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*BCE(i,j)*dy(i,j)-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                          
             elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) > 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                            
             elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+TE(i,j,n)*BCE(i,j)*dy(i,j)-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                          
             elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                          
             elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) == 100 && BCW(i,j) < 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                         
             elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) < 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                          
             elseif BCN(i,j) == 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) < 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) > 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                            
             elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) > 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) > 100
               
               P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                            
             elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) == 100
               
               P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                            
             elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) < 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*BCE(i,j)*dy(i,j)-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) < 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-TN(i,j,n)*BCN(i,j)*dx(i,j)+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                          
             elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) < 100
              
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*BCN(i,j)*dx(i,j)+TE(i,j,n)*BCE(i,j)*dy(i,j)-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                         
             elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) == 100
               
               P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-TN(i,j,n)*BCN(i,j)*dx(i,j)+TE(i,j,n)*BCE(i,j)*dy(i,j)-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                            
             elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) < 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) < 100
                
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                            
             elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) < 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+TE(i,j,n)*BCE(i,j)*dy(i,j)-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));                                                

                            
             elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+TE(i,j,n)*BCE(i,j)*dy(i,j)-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                            
             elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) > 100
                
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*BCE(i,j)*dy(i,j)-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                            
             elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) > 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-TN(i,j,n)*BCN(i,j)*dx(i,j)+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) > 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*BCN(i,j)*dx(i,j)+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                            
             elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-TN(i,j,n)*BCN(i,j)*dx(i,j)+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                          
             elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) < 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));

             elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) < 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                          
             elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) < 100
              
                P(i,j,n+1) = P(1i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) > 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*BCE(i,j)*dy(i,j)-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) > 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-TN(i,j,n)*BCN(i,j)*dx(i,j)+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                            
             elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) > 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*BCN(i,j)*dx(i,j)+TE(i,j,n)*BCE(i,j)*dy(i,j)-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-TN(i,j,n)*BCN(i,j)*dx(i,j)+TE(i,j,n)*BCE(i,j)*dy(i,j)-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                          
             elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+TE(i,j,n)*BCE(i,j)*dy(i,j)-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) == 100 && BCW(i,j) > 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+TE(i,j,n)*BCE(i,j)*dy(i,j)-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) > 100 && BCE(i,j) == 100 && BCS(i,j) < 100 && BCW(i,j) > 100
              
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-2*TN(i,j,n)*(P(i,j,n)-BCN(i,j))+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                           
             elseif BCN(i,j) == 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) > 100
                
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-2*TW(i,j,n)*(P(i,j,n)-BCW(i,j))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                            
             elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) == 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*BCS(i,j)*dx(i,j)-TN(i,j,n)*BCN(i,j)*dx(i,j)+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-TW(i,j,n)*(P(i,j,n)-P(i,j-1,n))+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));

             elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) == 100 && BCW(i,j) < 100
              
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(TS(i,j,n)*(P(i+1,j,n)-P(i,j,n))-TN(i,j,n)*BCN(i,j)*dx(i,j)+2*TE(i,j,n)*(BCE(i,j)-P(i,j,n))-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                            
             elseif BCN(i,j) < 100 && BCE(i,j) == 100 && BCS(i,j) > 100 && BCW(i,j) < 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-TN(i,j,n)*BCN(i,j)*dx(i,j)+TE(i,j,n)*(P(i,j+1,n)-P(i,j,n))-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
                          
             elseif BCN(i,j) == 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) < 100
               
                P(i,j,n+1) = P(i,j,n)+1/Cc(i,j)*(2*TS(i,j,n)*(BCS(i,j)-P(i,j,n))-TN(i,j,n)*(P(i,j,n)-P(i-1,j,n))+TE(i,j,n)*BCE(i,j)*dy(i,j)-TW(i,j,n)*BCW(i,j)*dy(i,j)+Qsc(i,j,n)-J(i,j)*BHP(i,j,n)*D(i,j)+Qg(i,j));
            
             elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block
               
             elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block 
 
             elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block                   

             elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
 
             elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block
 
             elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block
 
             elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block
 
             elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block
 
             elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block

             elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
 
             elseif BCN(i,j) < 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) > 100 %% Null block
 
             elseif BCN(i,j) > 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) < 100 %% Null block

             elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block 
 
             elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
  
             elseif BCN(i,j) < 100 && BCE(i,j) > 100 && BCS(i,j) > 100 && BCW(i,j) < 100 %% Null block
  
             elseif BCN(i,j) > 100 && BCE(i,j) < 100 && BCS(i,j) < 100 && BCW(i,j) > 100 %% Null block
 
             end
         end
     end
     
     if Cphi == 0
      
      break
      
  else
           
      if sum(sum(Pg1-P(:,:,n+1))) < Eps
          
          break
        
      else
      
      Pg1 = P(:,:,n+1);
      
      end
  end
  
  iter(n) = iteration;
  
 end
 
 if IBC == 1
 
 miu = ViscosityFunction(P(:,:,n+1),Nx,Ny,BCN,BCE,BCS,BCW); 
 Bo = FVFFunction(P(:,:,n+1),Nx,Ny,Co,B0,P0,BCN,BCE,BCS,BCW);
 J = (2*pi*1.127*(Kx.*Ky).^0.5.*dz./((log(0.2*dx./rw)+S).*miu(:,:,5).*Bo(:,:,5))); 
 BHP(:,:,n+1) = P(:,:,n+1)+Qsc(:,:,n)./J; 
 Qsc(:,:,n+1) = Qsc(:,:,n);
 Qqsc(:,:,n+1) = 0;
 
 for j = 1:Ny
            
     for i = 1:Nx
               
         if Qsc(i,j,n) == 0
                   
                    Qsc(i,j,n+1) = 0;
                    BHP(i,j,n+1) = 0;
                
         end         
     end     
 end
 
 else
     
        BHP(:,:,n+1) = BHP(:,:,n);
        Qqsc(:,:,n+1) = -J.*(P(:,:,n+1)-BHP(:,:,n+1));
        Qsc(:,:,n+1) = Qsc(:,:,n);
        
        for j = 1:Ny
        
            for i = 1:Nx
            
                if BHP(i,j,n) == 0
                
                    Qqsc(i,j,n+1)=0;
                    Qqsc(i,j,1)=0;
                
                end
            end
        end
 end
 
 MB = expMBCFunction(P(:,:,:),dx,dy,dz,dt,Nx,Ny,n,IBC,Qsc(:,:,:),Qqsc(:,:,:),Co,B0,P0,Phi0,Cphi,N_Aavg,E_Aavg,S_Aavg,W_Aavg,BCN,BCE,BCS,BCW);
 IMB(n) = MB(1);
 CMB(n) = MB(2);
 t(n+1) = sum(dt);

end

 
  

             
         
         
         
         
         
         
     
     
     
     
     
     
     
 
 
 
 
 
 
 
 
 

