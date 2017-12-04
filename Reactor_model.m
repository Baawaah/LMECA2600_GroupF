function [U5_burning_rate]= Reactor_model(t_final,n_th_init,mTot,U5_pour,U8_pour,Pu9_pour)
 close all, clear all
 
%% Données 
if nargin==0
    n_eV_th=0.025;
    n_eV_fast=10^6;
    Path='./DATABASE';
    X='U235';
    Transfo='Fission';
    t_final=10;
    n_th_init=10^10;
    mTot=25;
    U5_pour=3/100;
    U8_pour=97/100 ;
    Pu9_pour=0;
    V_core=10; %m3;
end
phi_n_th=(2200*n_th_init)/V_core;
phi_n_fast=0; 
% n_eV=0.025 %eV
bar_in=1; %pourcentage insertion des barres de contôle 
h=10^-2;  
dt=h;
z=0:h:t_final;
Avogadro=6.022*10^23;
%% Directory path 
Path='C:\Users\Pierre-Yves Legros\Documents\UCL\Nucléaire\Code\DATABASE';
%% Calcul des demi-vies et sections efficaces 



lambda_U9=log(2)/Demi_vie('U239','BetaMinus');
lambda_Np9=log(2)/Demi_vie('Np239','BetaMinus');
lambda_Pu9=log(2)/Demi_vie('Pu239','Alpha');
lambda_Pfp=log(2)/Demi_vie('Xe135','BetaMinus');
lambda_n=log(2)/(5*10^-4);
lambda_Pf=log(2);
lambda_rod = [50+bar_in*50,600+bar_in*1400] ;

sigma_U5f_th=10^-28*Section_efficace('U235','Fission',n_eV_th,Path);
sigma_U5f_fast=10^-28*Section_efficace('U235','Fission',n_eV_fast,Path);
sigma_U8f_th=10^-28*Section_efficace('U238','Fission',n_eV_th,Path);
sigma_U8f_fast=10^-28*Section_efficace('U238','Fission',n_eV_fast,Path);
sigma_U8a_th=10^-28*Section_efficace('U238','Capture',n_eV_th,Path);
sigma_U8a_fast=10^-28*Section_efficace('U238','Capture',n_eV_fast,Path);
sigma_U9f_th=10^-28*Section_efficace('U239','Fission',n_eV_th,Path);
sigma_U9f_fast=10^-28*Section_efficace('U239','Fission',n_eV_fast,Path);
sigma_Np9f_th=10^-28*Section_efficace('Np239','Fission',n_eV_th,Path);
sigma_Np9f_fast=10^-28*Section_efficace('Np239','Fission',n_eV_fast,Path);
sigma_Pu9f_th=10^-28*Section_efficace('Pu239','Fission',n_eV_th,Path);
sigma_Pu9f_fast=10^-28*Section_efficace('Pu239','Fission',n_eV_fast,Path);
sigma_Pfp_th=10^-28*Section_efficace('Xe135','Capture',n_eV_th,Path);
sigma_Pfp_fast=10^-28*Section_efficace('Xe135','Capture',n_eV_fast,Path);
coef_poison=0.05;



%% Fonction à résoudre
    function y =fun(~,x)
        phi_n_th=2200*x(6)/V_core ;
        phi_n_fast=1.4*10^7*x(7)/V_core;
        
        y=zeros(10,1);
        
        y(1)=lambda_Pu9*x(5)-(sigma_U5f_th*phi_n_th+sigma_U5f_fast*phi_n_fast)*x(1); %U5
        y(2)=-(sigma_U8f_th*phi_n_th+sigma_U8f_fast*phi_n_fast+sigma_U8a_th*phi_n_th+sigma_U8a_fast*phi_n_fast)*x(2); %U8
        y(3)=(sigma_U8a_th*phi_n_th+sigma_U8a_fast*phi_n_fast)*x(2)-(sigma_U9f_th*phi_n_th+sigma_U9f_fast*phi_n_fast)*x(3)-lambda_U9*x(3); %U9
        y(4)=-(sigma_Np9f_th*phi_n_th+sigma_Np9f_fast*phi_n_fast)*x(4)-lambda_Np9*x(4)+lambda_U9*x(3); %Np9
        y(5)=-(sigma_Pu9f_th*phi_n_th+sigma_Pu9f_fast*phi_n_fast)*x(5)+lambda_Np9*x(4)-lambda_Pu9*x(5); %Pu9
        y(6)=lambda_n*x(7)-(sigma_U5f_th*phi_n_th*x(1)+sigma_U8f_th*phi_n_th*x(2)+sigma_U9f_th*phi_n_th*x(3)+sigma_Np9f_th*phi_n_th*x(4)+sigma_Pu9f_th*phi_n_th*x(5))...
            -sigma_U8a_th*phi_n_th*x(2)+lambda_Pf*x(8)-lambda_rod(1)*x(6)-sigma_Pfp_th*phi_n_th*x(9); %n_th
        y(7)=-lambda_n*x(7)+(2*sigma_U5f_th*phi_n_th+sigma_U5f_fast*phi_n_fast)*x(1)+(2*sigma_U8f_th*phi_n_th+sigma_U8f_fast*phi_n_fast)*x(2)+(2*sigma_U9f_th*phi_n_th+sigma_U9f_fast*phi_n_fast)*x(3)...
            +(2*sigma_Np9f_th*phi_n_th+sigma_Np9f_fast*phi_n_fast)*x(4)+(2*sigma_Pu9f_th*phi_n_th+sigma_Pu9f_fast*phi_n_fast)*x(5)-sigma_U8a_fast*phi_n_fast*x(2)-lambda_rod(2)*x(7)-sigma_Pfp_fast*phi_n_fast*x(9);
        y(8)=(1-coef_poison)*((2*sigma_U5f_th*phi_n_th+2*sigma_U5f_fast*phi_n_fast)*x(1)+(2*sigma_U8f_th*phi_n_th+2*sigma_U8f_fast*phi_n_fast)*x(2)+(2*sigma_U9f_th*phi_n_th+2*sigma_U9f_fast*phi_n_fast)*x(3)...
            +(2*sigma_Np9f_th*phi_n_th+2*sigma_Np9f_fast*phi_n_fast)*x(4)+(2*sigma_Pu9f_th*phi_n_th+2*sigma_Pu9f_fast*phi_n_fast)*x(5))-lambda_Pf*x(8);
        y(9)=coef_poison*((2*sigma_U5f_th*phi_n_th+2*sigma_U5f_fast*phi_n_fast)*x(1)+(2*sigma_U8f_th*phi_n_th+2*sigma_U8f_fast*phi_n_fast)*x(2)+(2*sigma_U9f_th*phi_n_th+2*sigma_U9f_fast*phi_n_fast)*x(3)...
            +(2*sigma_Np9f_th*phi_n_th+2*sigma_Np9f_fast*phi_n_fast)*x(4)+(2*sigma_Pu9f_th*phi_n_th+2*sigma_Pu9f_fast*phi_n_fast)*x(5))-lambda_Pfp*x(9)-(sigma_Pfp_th*phi_n_th+sigma_Pfp_fast*phi_n_fast)*x(9);
        y(10)=lambda_Pf*x(8)+lambda_Pfp*x(9)+(sigma_Pfp_th*phi_n_th+sigma_Pfp_fast*phi_n_fast)*x(9);
       
    end



%% Conditions initiales

% Avogadro=1 ; 
 C_0=[(U5_pour*mTot)/(molarMass('U235'))*Avogadro,(U8_pour*mTot)/(molarMass('U238'))*Avogadro,0,0,(Pu9_pour*mTot)/(molarMass('Pu239'))*Avogadro,(10^10),0,0,0,0];
%% Méthode de résolution
 [t,C]=ode45(@(t,x) fun(t,x),z,C_0);
%   Avogadro=1 ; 
%% Plot des résultats 
% n=[0,t_final,10^-15,10^3]
loglog(t,C(:,5))


figure
loglog(t,C(:,1)/Avogadro)
hold on ;

loglog(t,C(:,2)/Avogadro)
hold on ;
loglog(t,C(:,3)/Avogadro)
hold on ;
loglog(t,C(:,4)/Avogadro)
hold on ;
loglog(t,C(:,5)/Avogadro)
hold on ;
loglog(t,C(:,6)/Avogadro)
hold on ;
loglog(t,C(:,7)/Avogadro)
hold on ;
loglog(t,C(:,8)/Avogadro)
hold on ;
loglog(t,C(:,9)/Avogadro)
hold on ;
loglog(t,C(:,10)/Avogadro)
% axis(n)
legend('U5','U8','U9','Np9','Pu9','n_{th}','n_{fast}','PF*','PF_p','PF')
hold off ; 



end