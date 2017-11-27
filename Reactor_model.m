function [U5_burning_rate]= Reactor_model(t_final,n_th_init,mTot,U5_pour,U8_pour,Pu9_pour)
 close all, clear all
if nargin==0
    n_eV_th=0.025;
    n_eV_fast=10^6;
    Path='./DATABASE';
    X='U235';
    Transfo='Fission';
    t_final=10;
    n_th_init=10^10;
    mTot=100;
    U5_pour=3/100;
    U8_pour=97/100 ;
    Pu9_pour=0;
    V_core=10; %m3;
end
phi_n_th=(2200*n_th_init)/V_core;
phi_n_fast=0; 
% n_eV=0.025 %eV

h=10^-2;  
dt=h;
z=0:h:t_final;
Avogadro=6.022*10^23;

Path='C:\Users\Pierre-Yves Legros\Documents\UCL\Nucléaire\Code\DATABASE';

lambda = log(2)./[Demi_vie('U239','BetaMinus'),Demi_vie('Np239','BetaMinus'),Demi_vie('Pu239','Alpha')];
sigma_th=10^-28.*[Section_efficace('U235','Fission',n_eV_th,Path),Section_efficace('U238','Fission',n_eV_th,Path),...
    Section_efficace('U238','Capture',n_eV_th,Path),Section_efficace('U239','Fission',n_eV_th,Path),Section_efficace('Np239','Fission',n_eV_th,Path)...
    Section_efficace('Pu239','Fission',n_eV_th,Path)];
sigma_fast = 10^-28*[Section_efficace('U235','Fission',n_eV_fast,Path),Section_efficace('U238','Fission',n_eV_fast,Path),...
    Section_efficace('U238','Capture',n_eV_fast,Path),Section_efficace('U239','Fission',n_eV_fast,Path),Section_efficace('Np239','Fission',n_eV_fast,Path)...
    Section_efficace('Pu239','Fission',n_eV_fast,Path)];


    function [y]= fun (~,x)
       phi_n_th=2200*x(7)/V_core ;
       phi_n_fast=1.4*10^7*x(8)/V_core;
       y=zeros(8,1);
       y(1)=-sigma_th(1)*phi_n_th*x(1)+lambda(3)*x(5)+2*(-sigma_fast(1)*phi_n_fast*x(1)+lambda(3)*x(5)); % U5
       y(2)=-sigma_th(2)*phi_n_th*x(2)-sigma_th(3)*phi_n_th*x(2)+2*(-sigma_fast(2)*phi_n_fast*x(2)-sigma_fast(3)*phi_n_fast*x(2)); %U8
       y(3)=sigma_th(3)*phi_n_th*x(2)-sigma_th(4)*phi_n_th*x(3)-lambda(1)*x(3)+2*(sigma_fast(3)*phi_n_fast*x(2)-sigma_fast(4)*phi_n_fast*x(3)-lambda(1)*x(3)); %U9
       y(4)=lambda(1)*x(3)-sigma_th(5)*phi_n_th*x(4)-lambda(2)*x(4)+2*(lambda(1)*x(3)-sigma_fast(5)*phi_n_fast*x(4)-lambda(2)*x(4)); %Np9
       y(5)=lambda(2)*x(4)-sigma_th(6)*phi_n_th*x(5)-lambda(3)*x(5)+2*(lambda(2)*x(4)-sigma_fast(6)*phi_n_fast*x(5)-lambda(3)*x(5)); %Pu9
       y(6)=sigma_th(1)*phi_n_th*x(1) + sigma_th(2)*phi_n_th*x(2) + sigma_th(4)*phi_n_th*x(3) + sigma_th(5)*phi_n_th*x(4) + sigma_th(6)*phi_n_th*x(5)...
           +2*(sigma_fast(1)*phi_n_fast*x(1) + sigma_fast(2)*phi_n_fast*x(2) + sigma_fast(4)*phi_n_fast*x(3) + sigma_fast(5)*phi_n_fast*x(4) + sigma_fast(6)*phi_n_fast*x(5)); %PF*
       y(7)=sigma_th(1)*phi_n_th*x(1)-sigma_th(3)*phi_n_th*x(2)+sigma_th(2)*phi_n_th*x(2)+sigma_th(4)*phi_n_th*x(3)+sigma_th(5)*phi_n_th*x(4)+sigma_th(6)*phi_n_th*x(5);  %n_th
       y(8)=sigma_fast(1)*phi_n_fast*x(1)+sigma_fast(2)*phi_n_fast*x(2)+sigma_fast(4)*phi_n_fast*x(3)+sigma_fast(5)*phi_n_fast*x(4)+sigma_fast(6)*phi_n_fast*x(5)%-sigma_fast(3)*phi_n_fast*x(2); %n_fast
    end

Avogadro=1;
 C_0=[(U5_pour*mTot)/(molarMass('U235'))*Avogadro,(U8_pour*mTot)/(molarMass('U238'))*Avogadro,0,0,(Pu9_pour*mTot)/(molarMass('Pu239'))*Avogadro,0,10^10,10^10];
 
 [t,C]=ode45(@(t,x) fun(t,x),z,C_0);
 
figure
subplot(3,3,1)
plot(t,C(:,1)/Avogadro)
title('Mole U5')
subplot(3,3,2)
plot(t,C(:,2)/Avogadro)
title('Mole U8')
subplot(3,3,3)
plot(t,C(:,3)/Avogadro)
title('Mole U9')
subplot(3,3,4)
plot(t,C(:,4)/Avogadro)
title('Mole Np9')
subplot(3,3,5)
plot(t,C(:,5)/Avogadro)
title('Mole Pu9')
subplot(3,3,6)
plot(t,C(:,6)/Avogadro)
title('Mole PF*')
subplot(3,3,7)
plot(t,C(:,7))
title('nombre neutron thermique')
subplot(3,3,8)
plot(t,C(:,8))
title('nombre neutron fast')
subplot(3,3,8)



end