% FLASH CALCULATION BY STABILITY ANALYSIS
clear
clc
close all

% Flash calculaTion for Twelve comPonenTs by sTabiliTy analize (PR) wiTh PR
% , SRK , SW



% +++++++++++++++++++++++++++++++  inPuTs  ++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

eos = input('1-PR 2-SRK 3-PT ?');
n=input('enTer number of comPonenTs =   ');
%for i=1:n   
 % z(i)=input('enter overal composition for component = ');
%end
% z= [0.3448 0.1379 0.05517 0.04138 0.06896 0.04827 0.02758 0.1379 0.006897 0.06896 0 0.6207];
 z=[0.3448 0.1379 0.0552 0.0414 0.0690 0.0483 0.0276 0.1379 0.0069 0.069 0.0021 0.06 ];
%P=input('enTer Pressure of sysTem in MPa =   ');
%T=input('enTer TemPreTure of sysTem in K =   ');
P=6.895;
T=344.3;
Tol=1E-12;
R=0.0083144;
delTa1=1+sqrt(2);
delTa2=1-sqrt(2);


kij = [0	0.0026	0.014	0.0256	0.0133	0.0056	0.0236	0.0422	0.0352	0.047	0.0474	0.05
0.0026	0	0.0011	0.0067	0.0096	0.008	0.0078	0.014	0.015	0.016	0.019	0.03
0.014	0.0011	0	0.0078	0.0033	0.0111	0.012	0.0267	0.056	0.059	0.007	0.02
0.0256	0.0067	0.0078	0	0	0.004	0.002	0.024	0.025	0.026	0.006	0.01
0.0133	0.0096	0.0033	0	0	0.017	0.017	0.0174	0.019	0.012	0.01	0.001
0.0056	0.008	0.0111	0.004	0.017	0	0	0	0	0	0	0
0.0236	0.0078	0.012	0.002	0.017	0	0	0	0	0	0	0
0.0422	0.014	0.0267	0.024	0.0174	0	0	0	0	0	0	0
0.0352	0.015	0.056	0.025	0.019	0	0	0	0	0	0	0
0.047	0.016	0.059	0.026	0.012	0	0	0	0	0	0	0
0.0474	0.019	0.007	0.006	0.01	0	0	0	0	0	0	0
0.05	0.03	0.02	0.01	0.001	0	0	0	0	0	0	0
];



Mw = xlsread('ProPerTy.xlsx',1,'B4:B15');
Tb = xlsread('ProPerTy.xlsx',1,'C4:C15');
Tc = xlsread('ProPerTy.xlsx',1,'D4:D15');
Pc = xlsread('ProPerTy.xlsx',1,'E4:E15');
Vc = xlsread('ProPerTy.xlsx',1,'F4:F15');
Zc = xlsread('ProPerTy.xlsx',1,'G4:G15');
acf = xlsread('ProPerTy.xlsx',1,'H4:H15');
Zra = xlsread('ProPerTy.xlsx',1,'I4:I15');
Sgr = xlsread('ProPerTy.xlsx',1,'J4:J15');
se = [-0.154 -0.1002 -0.08501 -0.07935 -0.06413 -0.0435 -0.04183 -0.01478 0 0 0 0];

for i=9:n
kwa(i)=4.5579*(Mw(i)^0.15178)*Sgr(i)^(-0.84573);
if kwa>12.5 & kwa<=13.5
    Se(i)= 1-2.258/(Mw(i)^0.1823);
elseif kwa>11 & kwa<=12.5
    Se(i)= 1-3.004/(Mw(i)^0.2324);
else  kwa>8.5 & kwa<=11;
      Se(i)= 1-2.516/(Mw(i)^0.2008);
    end
end

SE = Se + se;


% +++++++++++++++++++++++++++++++  PR ParameTers  +++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for i=1:n
    ac(i)=0.457235*(R^2)*((Tc(i))^2)/Pc(i);
    b(i)=0.077796*R*Tc(i)/Pc(i);
    m(i)=0.3796+1.485*acf(i)-0.1644*acf(i)^2+0.01677*acf(i)^3;
    alPha(i)=(1+m(i)*(1-sqrt(T/Tc(i))))^2;
    a(i)=alPha(i)*ac(i);
    c(i)=SE(i)*b(i);
end

% ++++++++++++++++++++++++  Mixing Rules  +++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

am=0;
for i=1:n
   amm = 0;
    for j=1:n
        amm=amm+z(i)*z(j)*sqrt(a(i)*a(j))*(1-kij(i,j));
    
    end
    am=amm+am;
end
bm = 0;
for i=1:n
bm = bm + z(i) * b(i);
end
A=am*P/((R*T)^2);
B=bm*P/(R*T);

% +++++++++++++++++++++++++++  Z EquaTion  ++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

U=2*B;
W=B;

zf=rooT_fcn(A,B,U,W);
zfl=min(zf);
zfv=max(zf);


% ++++++++++++++++++++++++  Psedou CriTical ProPerTies ++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% landa=(z'.*Vc)/(sum(z'.*Vc));
% Tcc=sum(landa.*Tc);
% TPc=sum(Tc.*z');
% PPc=sum(Pc.*z');
% wPc=sum(z'.*acf);
% Pcc=PPc*(1+(5.808+4.93*wPc)*((Tcc/TPc)-1));
% zPc=sum(z'.*Zc);
% Vcc=0.307*R*Tcc/Pcc;
% nucc=zPc*R*Tcc/Pcc;
Tpc=sum(Tc.*z');
er=1;
while er>10^(-10)
   TT=Tpc; 
for i=1:n
    ac(i)=0.457235*(R^2)*((Tc(i))^2)/Pc(i);
    b(i)=0.077796*R*Tc(i)/Pc(i);
    m(i)=0.3796+1.485*acf(i)-0.1644*acf(i)^2+0.01677*acf(i)^3;
    alPha(i)=(1+m(i)*(1-sqrt(TT/Tc(i))))^2;
    a(i)=alPha(i)*ac(i);
    c(i)=SE(i)*b(i);
end
am=0;
for i=1:n
   amm = 0;
    for j=1:n
        amm=amm+z(i)*z(j)*sqrt(a(i)*a(j))*(1-kij(i,j));
    
    end
    am=amm+am;
end
bm = 0;
for i=1:n
bm = bm + z(i) * b(i);
end
Tpc_new=(0.077796/0.457235)*(am/(R*bm));
er=abs(Tpc-Tpc_new);
Tpc=Tpc_new;
end
Ppc=0.077796*R*Tpc/bm;
Vcc=0.307*R*Tpc/Ppc;


% ++++++++++++++++++++++++++ K Wilson +++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for i=1:n
    if T>Tc(i)
   k(i)=(Pc(i)/P)*exp(5.37*(1+acf(i))*(1-Tc(i)/T))/(T/Tc(i));
    else
       k(i)=(Pc(i)/P)*exp(5.37*(1+acf(i))*(1-Tc(i)/T));  
    end

end

% +++++++++++++++++++++  fluid conditioning  ++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if zfl==zfv
    zf=zfv;
    Vcf=zf*R*T/P;
    
    if Vcf>Vcc
        like=2;
        Phase='gaslike';
        M=z./k;
        disp('gas like')
    elseif Vcf<Vcc
        like=1;
        Phase='liquidlike';
        M=z.*k;
        disp('liquid like')
    end
else        %% if zfl~=zfl
    %%gibbs
    %G=Gv-Gl
    
    G=(zfv-zfl)+log((zfl-B)/(zfv-B))-A/(B*(delTa2-delTa1))*log(((zfl+delTa1*B)*(zfv+delTa2*B))/((zfl+delTa2*B)*(zfv+delTa1*B)));
    if G<=0
        Phase='gaslike';
        disp('gas like')
        like=2;
        zf=zfv;
        M=z./k;
    else
        Phase='liquidlike';
        disp('liquid like')
        like=1;
        zf=zfl;
        M=z.*k;
    end
end



%% Phif
for i=1:n
    ac(i)=0.457235*(R^2)*((Tc(i))^2)/Pc(i);
    b(i)=0.077796*R*Tc(i)/Pc(i);
    m(i)=0.3796+1.485*acf(i)-0.1644*acf(i)^2+0.01677*acf(i)^3;
    alPha(i)=(1+m(i)*(1-sqrt(T/Tc(i))))^2;
    a(i)=alPha(i)*ac(i);
    c(i)=SE(i)*b(i);
end
am=0;
for i=1:n
    
    for j=1:n
        am=am+z(i)*z(j)*(((a(i)*a(j))^0.5)*(1-kij(i,j)));
    end
end
bm=sum(z.*b);
for i=1:n
    sumaij=0;
    for j=1:n
        sumaij=sumaij+(z(j)*(((a(i)*a(j))^0.5)*(1-kij(i,j))));
    end
     Phif(i)=exp((b(i)/bm)*(zf-1)-log(zf-B)-(A/(B*(delTa2-delTa1)))*...
        (((2*sumaij)/am)-(b(i)/bm))*log((zf+B*delTa2)/(zf+B*delTa1)));
    ff(i)=Phif(i)*z(i)*P;
end



%%%   sTabiliTy
e=1;
c=1;
while e>1e-12
    mi=M/sum(M);
    aM=0;
    for i=1:n
        for j=1:n
            aMij(i,j)=mi(i)*mi(j)*((a(i)*a(j))^0.5)*(1-kij(i,j));
            aM=aM+aMij(i,j);
        end
    end
    bM=sum(mi.*b);
    
    AM=aM*P/(R*T)^2;
    BM=bM*P/(R*T);
    UM=2*BM;
    WM=BM;

% zfM=rooT_fcn(AM,BM,UM,WM);
% zflM=min(zfM);
% zfvM=max(zfM);
%     
    if like==1    %% liquid like
        zM=max(rooT_fcn(AM,BM,UM,WM));
    elseif like==2    %% gas like
        zM=min(rooT_fcn(AM,BM,UM,WM));
    end
    
    
    for i=1:n
    ac(i)=0.457235*(R^2)*((Tc(i))^2)/Pc(i);
    b(i)=0.077796*R*Tc(i)/Pc(i);
    m(i)=0.3796+1.485*acf(i)-0.1644*acf(i)^2+0.01677*acf(i)^3;
    alPha(i)=(1+m(i)*(1-sqrt(T/Tc(i))))^2;
    a(i)=alPha(i)*ac(i);
    c(i)=SE(i)*b(i);
end
    
for i=1:n
    sumaij=0;
    for j=1:n
        sumaij=sumaij+(mi(j)*(((a(i)*a(j))^0.5)*(1-kij(i,j))));
    end
     PhiM(i)=exp((b(i)/bM)*(zM-1)-log(zM-BM)-(AM/(BM*(delTa2-delTa1)))*...
        (((2*sumaij)/aM)-(b(i)/bM))*log((zM+BM*delTa2)/(zM+BM*delTa1)));
%     ff(i)=PhiM(i)*z(i)*P;
end

fM=(PhiM).*(mi)*(P);

 %-log(PhiM(i))+log(zi(i)*Phif(i));%

% if like==1 % liquid like

k=Phif./PhiM;
M_new=(Phif./PhiM).*z;
e=sum((M_new-M).^2);

% 
% elseif like==2 % gas like
% %     k=fM./ff;
% %     e=sum((1-fM./ff).^2);
% %      MM=z./k;
% k=PhiM./Phif;
% M_new=(PhiM./Phif).*z;
% e=sum((M_new-M).^2);
% end

M=M_new;
c=c+1;

end
m=mi;
M;
sumM=sum(M);
% disP(['iniTial Phase is " ' Phase ' "'])
% 
% if sumM>1
%     disP('Feed is unsTable')
% else
%     disP('Feed is sTable')
% end
% 
if like ==1
    k=mi./z;
elseif like==2
    k=z./mi;
end

%%
if sumM > 1
error=1;

while error>10^(-12)
nv=rashfored(z,k);
yi=(z.*k)./(1+nv.*(k-1));
xi=z./(1+nv.*(k-1));

%%
for i=1:n
    ac(i)=0.457235*(R^2)*((Tc(i))^2)/Pc(i);
    b(i)=0.077796*R*Tc(i)/Pc(i);
    m(i)=0.3796+1.485*acf(i)-0.1644*acf(i)^2+0.01677*acf(i)^3;
    alPha(i)=(1+m(i)*(1-sqrt(T/Tc(i))))^2;
    a(i)=alPha(i)*ac(i);
    c(i)=SE(i)*b(i);
end
am=0;
for i=1:n
    
    for j=1:n
        am=am+xi(i)*xi(j)*(((a(i)*a(j))^0.5)*(1-kij(i,j)));
    end
end
bm=sum(xi.*b);

A=am*P/((R*T)^2);
B=bm*P/(R*T);

U=2*B;
W=B;

zf=rooT_fcn(A,B,U,W);
zf=min(zf);




for i=1:n
    sumaij=0;
    for j=1:n
        sumaij=sumaij+(xi(j)*(((a(i)*a(j))^0.5)*(1-kij(i,j))));
    end
     Phil(i)=exp((b(i)/bm)*(zf-1)-log(zf-B)-(A/(B*(delTa2-delTa1)))*...
        (((2*sumaij)/am)-(b(i)/bm))*log((zf+B*delTa2)/(zf+B*delTa1)));
    ffl(i)=Phil(i)*xi(i)*P;
end



%% 
% for i=1:n
%     ac(i)=0.457235*(R^2)*((Tc(i))^2)/Pc(i);
%     b(i)=0.077796*R*Tc(i)/Pc(i);
%     m(i)=0.3796+1.485*acf(i)-0.1644*acf(i)^2+0.01677*acf(i)^3;
%     alPha(i)=(1+m(i)*(1-sqrt(T/Tc(i))))^2;
%     a(i)=alPha(i)*ac(i);
%     c(i)=SE(i)*b(i);
% end

% ++++++++++++++++++++++++  Mixing Rules  +++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

for i=1:n
    ac(i)=0.457235*(R^2)*((Tc(i))^2)/Pc(i);
    b(i)=0.077796*R*Tc(i)/Pc(i);
    m(i)=0.3796+1.485*acf(i)-0.1644*acf(i)^2+0.01677*acf(i)^3;
    alPha(i)=(1+m(i)*(1-sqrt(T/Tc(i))))^2;
    a(i)=alPha(i)*ac(i);
    c(i)=SE(i)*b(i);
end
am=0;
for i=1:n
    
    for j=1:n
        am=am+yi(i)*yi(j)*(((a(i)*a(j))^0.5)*(1-kij(i,j)));
    end
end
bm=sum(yi.*b);
A=am*P/((R*T)^2);
B=bm*P/(R*T);

U=2*B;
W=B;

zf=rooT_fcn(A,B,U,W);
zf=max(zf);
for i=1:n
    sumaij=0;
    for j=1:n
        sumaij=sumaij+(yi(j)*(((a(i)*a(j))^0.5)*(1-kij(i,j))));
    end
     Phiv(i)=exp((b(i)/bm)*(zf-1)-log(zf-B)-(A/(B*(delTa2-delTa1)))*...
        (((2*sumaij)/am)-(b(i)/bm))*log((zf+B*delTa2)/(zf+B*delTa1)));
    ffv(i)=Phiv(i)*yi(i)*P;
end

knew=k.*(ffl./ffv);
k=knew;
error=sum((1-ffl./ffv).^2);


end
nl=1-nv;

elseif M < 1
    k;
    z;
end