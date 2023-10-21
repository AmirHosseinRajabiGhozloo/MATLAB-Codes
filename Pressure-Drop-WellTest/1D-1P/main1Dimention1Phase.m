clc;
clear;

%% Note : if RBC or LBC = 100 : boundary is located between two blocks and two blocks are interconnected
%% Note : if RBC or LBC < 100 : boundary is constant flow rate condition or dp/dx = c and the value of RBC or LBC is equal to c
%% Note : if RBC or LBC > 100 : boundary is constant pressure and at boundary we have P = C and the value of RBC or LBC is equal to C

%% ==========================================================================================================================================================
%% ++++++++++++++++++++++++++ INPUTS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ==========================================================================================================================================================
inputs = xlsread('Properties');
dx = inputs(:,1);                                                          %% the lenght of each grid
dy = inputs(:,2);                                                          %% the width of each grid fcvfd
dz = inputs(:,3);                                                          %% the height of each grid
Ngrid = numel(dx);                                                         %% How many grids do we have(need) ?
LBC = inputs(:,4);                                                         %% Left Boundary Condition of each grid
RBC = inputs(:,5);                                                         %% Right Boundary Condition of each grid
initial_BC = inputs(1,6);                                                  %% Constant Flow Rate or Constant BHP ?
BHP = inputs(:,7);                                                         %% Bottem Hole Pressure
Qsc = inputs(:,8);                                                         %% well Flow Rates
rw = inputs(:,9);                                                          %% Well Radiouses
S = inputs(:,10);                                                          %% Skin Factor
dt = inputs(1,11);                                                         %% each Time Steps Duration
Tt = inputs(1,12);                                                         %% Total time used for Simulation
t(1) = 0;                                                                  %% starting time
Np(1) = 0;                                                                 %% total production at the starting time
Pi = inputs(1,13);                                                         %% Initial Pressure
Co = inputs(1,14);                                                         %% Oil Comperssibility
B0 = inputs(1,15);                                                         %% Refrence Formation Volume Factor
Niteration = input('How many Iterations Do you Prefer to Run ?');          %% Number of Iterations Run
Solvation = inputs(1,16);                                                  %% Implicit or Explicit ?
if Solvation == 1
    Ntime = (Tt/dt);                                                       %% number of time steps for Implicit Solver
else
    Ntime = input('How many Time Steps Do You Prefer ?');                  %% number of time steps for Explicit Solver
end
Eps = input('How much Pressure diffrence is acceptable for breaking loop ?'); 
Rho = inputs(1,17);                                                        %% Oil Density
Gama = (0.00021584*(Rho)*32.17);                                           %% specific Density
Top = inputs(:,18);                                                        %% the Depth of Top of each grid
Zmid = Top + (dz/2);                                                       %% the Depth of middle of each grid
Kx = inputs(:,19);                                                         %% Permeability of each grid in x direction
Phi0 = inputs(1,20);                                                       %% Refrence porosity
Cphi = inputs(1,21);                                                       %% Rock Comperssibility
P(1:Ngrid,1) = Pi;                                                         %% the matrix of Pressure
P0 = inputs(1,22);                                                         %% Refrence Pressure
Qsc(:,1) = Qsc;                                                            %% the matrix of Flow Rate
Aavg = AavgFunction(dx,dy,dz,Ngrid,Kx,LBC,RBC);                            %% Aavg Function
L_Aavg = Aavg(:,1);                                                        %% left Average flow Area for each grid
R_Aavg = Aavg(:,2);                                                        %% Right Average flow Area for each grid
req = (((1/pi)^(0.5))*(((dx.^2)+(dy.^2)).^(0.5)));                         %% equal Radious for each grid
x(1) = ((dx(1))/2);                                                          

for i=1:Ngrid-1                                                            %% x-Location of middle of each grid
    
    x(i+1) = (x(i)+((dx(i)+dx(i+1))/2));
    
end

K = PermFunction(x,Ngrid,Kx);                                              %% Permeability based on location
miu = ViscosityFunction(P(:,1),Ngrid,LBC,RBC);                             %% Viscosity of each grid
Boi = FVFFunction(Ngrid,P(:,1),Co,B0,P0,LBC,RBC);                          %% Formation Voleum Factor for each grid 
J = (2*pi*1.127.*Kx.*dz./(miu(:,3).*Boi(:,3).*(log(req./rw)+S)));          %% Productivity index

if initial_BC == 2
    
    Qqsc(:,1) = (-J.*(P(:,1)-BHP(1,:)'));
    
end

%% ==========================================================================================================================================================
%% +++++++++++++++++ FUNCTION'S UPDATE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ==========================================================================================================================================================

for j = 1:Ntime
    
    miu = ViscosityFunction(P(:,j),Ngrid,LBC,RBC);
    Lmiu(:,j) = miu(:,1);                                                  %% viscosity at Left Boundary Condition
    Rmiu(:,j) = miu(:,2);                                                  %% viscosity at Right Boundary Condition                                                
    Bo = FVFFunction(Ngrid,P(:,j),Co,B0,P0,LBC,RBC);
    LBo(:,j) = Bo(:,1);                                                    %% Formation Volume Factor at Left Boundary Condition 
    RBo(:,j) = Bo(:,2);                                                    %% Formation Volume Factor at Right Boundary Condition 
    Phi(:,j) = PorosityFunction(P(:,j),Ngrid,P0,Phi0,Cphi);                     %% Porosity at each grid
    LTx(:,j) = 1.127*L_Aavg./(Lmiu(:,j).*LBo(:,j));                        %% Transmissibility at Left boundary
    RTx(:,j) = 1.127*R_Aavg./(Rmiu(:,j).*RBo(:,j));                        %% Transmissibility at Right boundary
    Qg = GravityFunction(Ngrid,Zmid,Gama,LBC,RBC,LTx,RTx);                 %% Flow because of grid Elevation
    J = (2*pi*1.127.*Kx.*dz./(miu(:,3).*Boi(:,3).*(log(req./rw)+S)));      %% Productivity index
    Pg = P(:,j);                                                           %% Initial Prssure Guess For Iterations
    
    for Iteration = 1:Niteration
        
     Boi = FVFFunction(Ngrid,Pg,Co,B0,P0,LBC,RBC);  
     
%% ==========================================================================================================================================================
%% +++++++++++++++ SOLVATION = 1  &&  INITIAL_BC = 1 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ==========================================================================================================================================================
     
     if Solvation == 1
         
         if initial_BC == 1
             
           Cc = dx.*dy.*dz/(5.615*dt).*(Phi(:,j)*Co/B0+Cphi*Phi0./Boi(:,3));  %% The Factor of Right side of the Equation
           
           for i = 1:Ngrid
           
           if LBC(i) == 100 && RBC(i) == 100
               
                 A(i,i,j) = (-(RTx(i,j)+LTx(i,j))-Cc(i));
                 A(i,i+1,j) = RTx(i,j);
                 A(i,i-1,j) = LTx(i,j);
                 B(i,j) = (-Qsc(i,j)-Cc(i)*P(i,j)+Qg(i));
             
           elseif LBC(i) == 100 && RBC(i) < 100
               
                 A(i,i,j) = (-LTx(i,j)-Cc(i));
                 A(i,i-1,j) = LTx(i,j);
                 B(i,j) = (-RTx(i)*RBC(i)*dx(i)-Qsc(i,j)-Cc(i)*P(i,j)+Qg(i));
                 
           elseif LBC(i) == 100 && RBC(i) > 100   
               
                 A(i,i,j) = (-(2*RTx(i,j)+LTx(i,j))-Cc(i));
                 A(i,i-1,j) = LTx(i,j);
                 B(i,j) = (-2*RTx(i,j)*RBC(i)-Qsc(i,j)-Cc(i)*P(i,j)+Qg(i));
                 
           elseif LBC(i) < 100 && RBC(i) == 100
               
                 A(i,i,j) = (-RTx(i,j)-Cc(i));
                 A(i,i+1,j) = RTx(i,j);
                 B(i,j) = (LTx(i,j)*LBC(i)*dx(i)-Qsc(i,j)-Cc(i)*P(i,j)+Qg(i));
                 
           elseif LBC(i) > 100 && RBC(i) == 100
               
                 A(i,i,j) = (-(RTx(i,j)+2*LTx(i,j))-Cc(i));
                 A(i,i+1,j) = RTx(i,j);
                 B(i,j) = (-2*LTx(i,j)*LBC(i)-Qsc(i,j)-Cc(i)*P(i,j)+Qg(i));
                 
           end
           end
           
           
            P(:,j+1) = inv(A(:,:,j))*B(:,j);
            
            for i = 1:Ngrid
                
                Qsc(i,j+1) = Qsc(i,j);
                Qqsc(i) = 0;
                
                if Qsc(i,j) == 0
                    
                    BHP(i,j+1) = 0;
                    
                else 
                    
                    BHP(i,j+1) = Qsc(i,j)/J(i)+P(i,j+1);
                    
                end
            end
            
%% ==========================================================================================================================================================
%% ++++++++++++++++ SOLVATION = 1  &&  INITIAL_BC = 2 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ==========================================================================================================================================================
    
         elseif initial_BC == 2
             
             Cc = dx.*dy.*dz./(5.615*dt).*(Phi(:,j)*Co/B0+Cphi*Phi0./Boi(:,3));
             
             for i = 1:Ngrid
                 
                 if LBC(i) == 100 && RBC(i) == 100
                     
                     if BHP(i) == 0
                         
                          A(i,i,j) = (-(RTx(i,j)+LTx(i,j))-Cc(i));
                          A(i,i+1,j) = RTx(i,j);
                          A(i,i-1,j) = LTx(i,j);
                          B(i,j) = (-Cc(i)*P(i,j)+Qg(i));
                          
                     else 
                         
                          A(i,i,j) = (-(RTx(i,j)+LTx(i,j)+J(i))-Cc(i));
                          A(i,i+1,j) = RTx(i,j);
                          A(i,i-1,j) = LTx(i,j);
                          B(i,j) = (-Cc(i)*P(i,j)-J(i)*BHP(i,j)+Qg(i));
                          
                     end
                     
                 elseif LBC(i) == 100 && RBC(i) < 100
                     
                     if BHP(i) == 0
                         
                          A(i,i,j) = (-LTx(i,j)-Cc(i));
                          A(i,i-1,j) = LTx(i,j);
                          B(i,j) = (-RTx(i)*RBC(i)*dx(i)-Cc(i)*P(i,j)+Qg(i));
                          
                     else
                         
                          A(i,i,j) = (-LTx(i,j)-J(i)-Cc(i));
                          A(i,i-1,j) = LTx(i,j);
                          B(i,j) = (-RTx(i)*RBC(i)*dx(i)-J(i)*BHP(i,j)-Cc(i)*P(i,j)+Qg(i));
                          
                     end
                     
                 elseif LBC(i) == 100 && RBC(i) > 100
                     
                     if BHP(i) == 0
                         
                          A(i,i,j) = (-(2*RTx(i,j)+LTx(i,j))-Cc(i));
                          A(i,i-1,j) = LTx(i,j);
                          B(i,j) = (-2*RTx(i,j)*RBC(i)-Cc(i)*P(i,j)+Qg(i));
                          
                     else
                         
                          A(i,i,j) = (-(2*RTx(i,j)+LTx(i,j)+J(i))-Cc(i));
                          A(i,i-1,j) = LTx(i,j);
                          B(i,j) = (-2*RTx(i,j)*RBC(i)-Qsc(i,j)-J(i)*BHP(i,j)-Cc(i)*P(i,j)+Qg(i));
                          
                     end
                     
                 elseif LBC(i) < 100 && RBC(i) == 100
                     
                     if BHP(i) == 0
                         
                          A(i,i,j) = (-RTx(i,j)-Cc(i));
                          A(i,i+1,j) = RTx(i,j);
                          B(i,j) = (LTx(i,j)*LBC(i)*dx(i)-Cc(i)*P(i,j)+Qg(i));
                          
                     else
                         
                          A(i,i,j) = (-RTx(i,j)-J(i)-Cc(i));
                          A(i,i+1,j) = RTx(i,j);
                          B(i,j) = (LTx(i,j)*LBC(i)*dx(i)-J(i)*BHP(i,j)-Cc(i)*P(i,j)+Qg(i));
                         
                     end
                     
                 elseif LBC(i) > 100 && RBC(i) == 100
                     
                     if BHP(i) == 0
                         
                         A(i,i,j) = (-(RTx(i,j)+2*LTx(i,j))-Cc(i));
                         A(i,i+1,j) = RTx(i,j);
                         B(i,j) = (-2*LTx(i,j)*LBC(i)-Cc(i)*P(i,j)+Qg(i));
                         
                     else
                         
                         A(i,i,j) = (-(RTx(i,j)+2*LTx(i,j)+J(i))-Cc(i));
                         A(i,i+1,j) = RTx(i,j);
                         B(i,j) = (-2*LTx(i,j)*LBC(i)-J(i)*BHP(i,j)-Cc(i)*P(i,j)+Qg(i));
                         
                     end
                 end
             end
             
             P(:,j+1) = inv(A(:,:,j))*B(:,j);
             Qqsc(:,j+1) = (-J.*(P(:,j+1)-BHP(:,j)));
             
             for i = 1:Ngrid
             
                 if BHP(i,j) == 0
                 
                     Qqsc(i,j+1)=0;
                
                 end
                 
                 BHP(i,j+1)=BHP(i,j);
                
             end
             
             Np(j+1)=Np(j)+sum(abs(dt*Qqsc(:,j+1)));
             RF(j+1)=Np(j+1)*5.615/sum(Boi(:,3).*dx.*dy.*dz.*Phi(:,1));
         
         end
         
%% ==========================================================================================================================================================
%% +++++++++++++ SOLVATION = 2  &&  INITIAL_BC = 1 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ==========================================================================================================================================================
             
     elseif Solvation == 2
         
         if initial_BC == 1
             
             Landa = (dx.*dy.*dz./(5.615).*(Phi(:,j)*Co/B0+Cphi*Phi0./Boi(:,3)));
             F = Landa./(RTx(:,j)+LTx(:,j));
             dt(j) = min(F);
             Cc = dx.*dy.*dz./(5.615*dt(j)).*(Phi(:,j)*Co/B0+Cphi*Phi0./Boi(:,3));
             
             for i = 1:Ngrid
                 
                 if LBC(i) == 100 && RBC(i) == 100
                     
                     P(i,j+1) = (P(i,j)+1/Cc(i)*Qsc(i,j)+1/Cc(i)*(RTx(i,j)*P(i+1,j)-(RTx(i,j)+LTx(i,j))*P(i,j)+LTx(i,j)*P(i-1,j)-Qg(i)));
                      
                 elseif LBC(i) == 100 && RBC(i) < 100
                     
                     P(i,j+1) = (P(i,j)+1/Cc(i)*Qsc(i,j)+1/Cc(i)*(RTx(i)*RBC(i)*dx(i)-LTx(i,j)*P(i,j)+LTx(i,j)*P(i-1,j)-Qg(i)));
                     
                 elseif LBC(i) == 100 && RBC(i) > 100
                     
                     P(i,j+1) = (P(i,j)+1/Cc(i)*Qsc(i,j)+1/Cc(i)*(2*RTx(i,j)*RBC(i)-(2*RTx(i,j)+LTx(i,j))*P(i,j)+LTx(i,j)*P(i-1,j)-Qg(i)));  
                     
                 elseif LBC(i) < 100 && RBC(i) == 100
                     
                      P(i,j+1) = (P(i,j)+1/Cc(i)*Qsc(i,j)+1/Cc(i)*(RTx(i,j)*P(i+1,j)-RTx(i,j)*P(i,j)-LTx(i)*LBC(i)*dx(i)-Qg(i)));
                      
                 elseif LBC(i) > 100 && RBC(i) == 100
                     
                     P(i,j+1) = (P(i,j)+1/Cc(i)*Qsc(i,j)+1/Cc(i)*(RTx(i,j)*P(i+1,j)-(RTx(i,j)+2*LTx(i,j))*P(i,j)+2*LTx(i,j)*LBC(i)-Qg(i)));
                     
                 end
             end
             
              for i = 1:Ngrid
            
                  Qsc(i,j+1) = Qsc(i,j);
                  Qqsc(i) = 0;
                  
                  if Qsc(i,j) == 0
                 
                      BHP(i,j+1) = 0;
             
                  else
                      
                      BHP(i,j+1) = Qsc(i)/J(i)+P(i,j+1);
             
                  end
              end
              
%% ==========================================================================================================================================================
%% ++++++++++ SOLVATION = 2  &&  INITIAL_BC = 2 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ==========================================================================================================================================================
              
         elseif initial_BC == 2
             
             Landa(:,j) = (dx.*dy.*dz./(5.615).*(Phi(:,j)*Co/B0+Cphi*Phi0./Boi(:,3)));
             F(:,j) = Landa(:,j)./(RTx(:,j)+LTx(:,j)+J);
             dt(j) = min(F(:,j));
             Cc = (dx.*dy.*dz./(5.615*dt(j)).*(Phi(:,j)*Co/B0+Cphi*Phi0./Boi(:,3)));
             
             for i = 1:Ngrid
                 
                 if LBC(i) == 100 && RBC(i) == 100
                     
                     if BHP(i) == 0
                         
                       P(i,j+1) = (P(i,j)+1/Cc(i)*(RTx(i,j)*P(i+1,j)-(RTx(i,j)+LTx(i,j))*P(i,j)+LTx(i,j)*P(i-1,j)-Qg(i))); 
                       
                     else
                         
                       P(i,j+1) = (P(i,j)-1/Cc(i)*J(i)*(P(i,j)-BHP(i,j))+1/Cc(i)*(RTx(i,j)*P(i+1,j)-(RTx(i,j)+LTx(i,j))*P(i,j)+LTx(i,j)*P(i-1,j)-Qg(i)));
                         
                     end
                     
                 elseif LBC(i) == 100 && RBC(i) < 100
                     
                     if BHP(i) == 0
                         
                       P(i,j+1) = (P(i,j)+1/Cc(i)*(RTx(i)*RBC(i)*dx(i)-LTx(i,j)*P(i,j)+LTx(i,j)*P(i-1,j)-Qg(i)));   
                         
                     else
                         
                       P(i,j+1) = (P(i,j)-1/Cc(i)*J(i)*(P(i,j)-BHP(i,j))+1/Cc(i)*(RTx(i)*RBC(i)*dx(i)-LTx(i,j)*P(i,j)+LTx(i,j)*P(i-1,j)-Qg(i)));
                         
                     end
                     
                 elseif LBC(i) == 100 && RBC(i) > 100
                     
                     if BHP(i) == 0
                         
                       P(i,j+1) = (P(i,j)+1/Cc(i)*(2*RTx(i,j)*RBC(i)-(2*RTx(i,j)+LTx(i,j))*P(i,j)+LTx(i,j)*P(i-1,j)-Qg(i)));  
                   
                     else
                         
                       P(i,j+1) = (P(i,j)-1/Cc(i)*J(i)*(P(i,j)-BHP(i,j))+1/Cc(i)*(2*RTx(i,j)*RBC(i)-(2*RTx(i,j)+LTx(i,j))*P(i,j)+LTx(i,j)*P(i-1,j)-Qg(i)));
                     
                     end
                     
                 elseif LBC(i) < 100 && RBC(i) == 100
                     
                     if BHP(i) == 0
                         
                       P(i,j+1) = (P(i,j)+1/Cc(i)*(RTx(i,j)*P(i+1,j)-RTx(i,j)*P(i,j)-LTx(i)*LBC(i)*dx(i)-Qg(i)));
                  
                     else
                         
                       P(i,j+1) = (P(i,j)-1/Cc(i)*J(i)*(P(i,j)-BHP(i,j))+1/Cc(i)*(RTx(i,j)*P(i+1,j)-RTx(i,j)*P(i,j)-LTx(i)*LBC(i)*dx(i)-Qg(i)));
                         
                     end
                     
                 elseif LBC(i) > 100 && RBC(i) == 100
                     
                     if BHP(i) == 0
                         
                        P(i,j+1) = (P(i,j)+1/Cc(i)*(RTx(i,j)*P(i+1,j)-(RTx(i,j)+2*LTx(i,j))*P(i,j)+2*LTx(i,j)*LBC(i)-Qg(i))); 
                     
                     else
                         
                        P(i,j+1) = (P(i,j)-1/Cc(i)*J(i)*(P(i,j)-BHP(i,j))+1/Cc(i)*(RTx(i,j)*P(i+1,j)-(RTx(i,j)+2*LTx(i,j))*P(i,j)+2*LTx(i,j)*LBC(i)-Qg(i)));
                
                     end
                 end
             end
             
          Qqsc(:,j+1) = (-J.*(P(:,j+1)-BHP(:,j)));
                
          for i = 1:Ngrid
                
              if BHP(i,j) == 0
                    
                  Qqsc(i,j+1) = 0;
               
              end
              
              BHP(i,j+1) = BHP(i,j);
               
          end
          
          Np(j+1) = Np(j)+sum(abs(dt(j)*Qqsc(:,j+1)));
          RF(j+1) = Np(j+1)*5.615/sum(Boi(:,3).*dx.*dy.*dz.*Phi(:,1));
          
         end
     end
     
%% ==========================================================================================================================================================
%% ++++++++++++ PRESSURE AND TIME UPDATE ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ==========================================================================================================================================================
    
      if Solvation == 2
        
         t(j+1) = t(j)+dt(j);
       
     else
         
         t(j+1) = t(j)+dt;
       
     end
     
     if Cphi == 0
            
         break
       
     else
         
         if sum(abs(P(:,j+1)-Pg)) < Eps
            
             break
        
         else
             
             Pg = P(:,j+1);
       
         end
     end
    end
    
%% ==========================================================================================================================================================
%% +++++++++++++++ MATERIAL BALANCE CONTROL +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ==========================================================================================================================================================
  
    if Solvation == 2
        
       MB = MBCFunction(dx,dy,dz,dt(j),Ngrid,LBC,RBC,L_Aavg,R_Aavg,P(:,:),P0,Co,B0,Phi0,Cphi,j,Qsc(:,:),Qqsc(:,:),initial_BC); 

   else
        
       MB = MBCFunction(dx,dy,dz,dt,Ngrid,LBC,RBC,L_Aavg,R_Aavg,P(:,:),P0,Co,B0,Phi0,Cphi,j,Qsc(:,:),Qqsc(:,:),initial_BC); 
    
    end

   IMB(j) = MB(1);
   CMB(j) = MB(2);
 
end

%% ==========================================================================================================================================================
%% +++++++++++++++++ THE END ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% ==========================================================================================================================================================

