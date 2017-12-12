function [Lambda_BC_thermal,Lambda_BC_fast,U5_burning_rate]= Reactor_model(t_final,n_th_init,mTot,U5_pour,U8_pour,Pu9_pour,Poisson_pourc,mode)
 close all
 
 global  Power_old bar_in;
%% Données 
if nargin==0
    n_eV_th=0.025;
    n_eV_fast=10^6;
    Path='./DATABASE';
    X='U235';
    Transfo='Fission';
    t_final=500;
    n_th_init=10^10;
    mTot=25;
    U5_pour=3/100;
    U8_pour=97/100 ;
    Pu9_pour=0;
    V_core=10; %m3;
    Poisson_pourc= 0.05 ; 
    mode=3;
end
bar_in=0,1896 ; 
bar_in_pal=0.1896 ;
phi_n_th=(2200*n_th_init)/V_core;
phi_n_fast=0;   
dt=10^-4;
z=0:dt:t_final;
Power_target=10^4;
Power_cible=linspace(1000,10^5,3)
Time=linspace(100,t_final-100,3)
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
coef_poison=Poisson_pourc;

save=0;
Power_old=0;
k=0;


    function [P]= findpower(t)
        for i=1:(length(Time)-1)
            if t>=Time(i) && t<Time(i+1)
                P=Power_cible(i);
                return
            end
        end
        P=Power_cible(end);
        
    end


%% Fonction à résoudre
    function [y] =fun(t,x)
        
        lambda_rod =[50+bar_in*50,600+bar_in*1400];
        
        phi_n_th=2200*x(6)/V_core ;
        phi_n_fast=1.4*10^7*x(7)/V_core;
        
        y=zeros(11,1);
        
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
            +(2*sigma_Np9f_th*phi_n_th+2*sigma_Np9f_fast*phi_n_fast)*x(4)+(2*sigma_Pu9f_th*phi_n_th+2*sigma_Pu9f_fast*phi_n_fast)*x(5))-lambda_Pf*x(8); %PF*
        y(9)=coef_poison*((2*sigma_U5f_th*phi_n_th+2*sigma_U5f_fast*phi_n_fast)*x(1)+(2*sigma_U8f_th*phi_n_th+2*sigma_U8f_fast*phi_n_fast)*x(2)+(2*sigma_U9f_th*phi_n_th+2*sigma_U9f_fast*phi_n_fast)*x(3)...
            +(2*sigma_Np9f_th*phi_n_th+2*sigma_Np9f_fast*phi_n_fast)*x(4)+(2*sigma_Pu9f_th*phi_n_th+2*sigma_Pu9f_fast*phi_n_fast)*x(5))-lambda_Pfp*x(9)-(sigma_Pfp_th*phi_n_th+sigma_Pfp_fast*phi_n_fast)*x(9); %PF*_p
        y(10)=lambda_Pf*x(8)+lambda_Pfp*x(9)+(sigma_Pfp_th*phi_n_th+sigma_Pfp_fast*phi_n_fast)*x(9);
        
        Power_now =((x(8)+x(9))*200*10^6+lambda_Pf*x(8)*5*10^6+0.025*10^-6*x(7)*lambda_n)/dt*1.60218e-19;
        if mode == 1

            error=Power_now-Power_old;
            k=0.05*dt ; %5 max par seconde 
            if abs(error) > 10^-12 && error > 0
                bar_in = min(bar_in+k,1) ; 
            elseif abs(error) > 10^-12 && error < 0
                bar_in=max(bar_in-k,0);
            end
            Power_old = Power_now ;

            y(11)=bar_in ; 
        elseif mode ==2
        
          
            error=Power_now-Power_target ;
            k=0.05*dt ; %5 max par seconde 
            if abs(error) > 10^-12 && error > 0
                bar_in = min(bar_in+k,1) ; 
            elseif abs(error) > 10^-12 && error < 0
                bar_in=max(bar_in-k,0);
            end
            Power_old = Power_now ;

            y(11)=bar_in ; 

            
        elseif mode ==3 
            
            error=findpower(t)-Power_now ;
            k=0.05*dt; %5 max par seconde 
            if abs(error) > 10^-4 && error > 0
                bar_in = min(bar_in+k,0.6) ; 
            elseif abs(error) > 10^-4 && error < 0
                bar_in=max(bar_in-k,0);
            end
            Power_old = Power_now ;

            y(11)=bar_in ; 
        end
        
    end

%% Conditions initiales

% Avogadro=1 ; 
 C_0=[(U5_pour*mTot)/(molarMass('U235'))*Avogadro,(U8_pour*mTot)/(molarMass('U238'))*Avogadro,0,0,(Pu9_pour*mTot)/(molarMass('Pu239'))*Avogadro,(10^10),0,0,0,0,0];
%% Méthode de résolution
 [t,C]=ode45(@(t,x) fun(t,x),z,C_0);
 Power_reactor = ((C(:,8)+C(:,9))*200*10^6+lambda_Pf*C(:,8)*5*10^6+0.025*10^-6*C(:,7)*lambda_n)/dt*1.60218e-19;
 Power_split = [(C(:,8)+C(:,9))*200*10^6/dt,lambda_Pf*C(:,8)*5*10^6/dt,0.025*10^-6*C(:,7)*lambda_n/dt]*1.60218e-19;
%   Avogadro=1 ; 

%% Calcul des outputs

a=(diff((C(:,11)))/dt);
bar_critic=a(end)
Lambda_BC_thermal=50+bar_critic*50;
Lambda_BC_fast=600+bar_critic*1400;

U5_burning_rate=(C(1,1)-C(end,1))/C(1,1);

%% Plot des résultats 
n=[100,t_final,-Inf,Inf]



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

figure
loglog(t,Power_reactor)
axis(n)
hold on ; 
figure
pie (Power_split(end,:))
legend('Fission','Stabilization','Neutron conversion')
hold on;
end