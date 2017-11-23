function [U5_burning_rate]= Reactor_model(t_final,n_th_init,mTot,U5_pour,U8_pour,Pu9_pour)

if nargin==0
    n_eV=logspace(-5,6,10000);
    Path='./DATABASE';
    X='U235';
    Transfo='Fission';
    t_final=10^-2;
    n_th_init=10^10;
    mTot=100;
    U5_pour=3/100;
    U8_pour=97/100 ;
    Pu9_pour=0;
    
end
phi_n=2200*n_th_init*10^-4;
n_eV=0.025 %eV

h=10^-2;   
% z=0:h:t_final;
z=linspace(0,10^-2,5)


lambda = log(2)./[Demi_vie('U239','BetaMinus'),Demi_vie('Np239','BetaMinus'),Demi_vie('Pu239','Alpha')];
sigma=10^-24.*[Section_efficace('U235','Fission',n_eV,Path),Section_efficace('U238','Fission',n_eV,Path),...
    Section_efficace('U238','Capture',n_eV,Path),Section_efficace('U239','Fission',n_eV,Path),Section_efficace('Np239','Fission',n_eV,Path)...
    Section_efficace('Pu239','Fission',n_eV,Path)];

sigma
Path='C:\Users\Pierre-Yves Legros\Documents\UCL\Nucléaire\Code\DATABASE';

    function [y]= fun (t,x)
        
        y=zeros(5,1);
       y(1)=-sigma(1)*phi_n*x(1)+lambda(3)*x(5);
       y(2)=-sigma(2)*phi_n*x(2)-sigma(3)*phi_n*x(2);
       y(3)=sigma(3)*phi_n*x(2)-sigma(4)*phi_n*x(3)-lambda(1)*x(3);
       y(4)=lambda(1)*x(3)-sigma(5)*phi_n*x(4)-lambda(2)*x(4);
       y(5)=lambda(2)*x(4)-sigma(6)*phi_n*x(5)-lambda(3)*x(5);
       y(6)=-sigma(3)*phi_n*x(2);
       y(7)=sigma(1)*phi_n*x(1) + sigma(2)*phi_n*x(2) + sigma(4)*phi_n*x(3) + sigma(5)*phi_n*x(4) + sigma(6)*phi_n*x(5);
       
       %phi_n = 2200*(x(6) + y(6));
       y
    end

 C_0=[(U5_pour*mTot)/molarMass('U235'),(U8_pour*mTot)/molarMass('U238'),0,0,(Pu9_pour*mTot)/molarMass('Pu239'),n_th_init,0];
X=ode45(@fun,z,C_0)
plot(X.x,X.y)
end


%%
% % It calculates ODE using Runge-Kutta 4th order method
% y=zeros(length(z),5);
% C_0=[(U5_pour*mTot)/molarMass('U235'),(U8_pour*mTot)/molarMass('U238'),0,0,(Pu9_pour*mTot)/molarMass('Pu239')];
% y(1,:)=C_0;
% 
% for i=1:(length(z)-1)                              % calculation loop
%     k_1 = fun(z(i),y(i,:));
%     k_2 = fun(z(i)+0.5*h,y(i,:)+0.5*h*k_1);
%     k_3 = fun((z(i)+0.5*h),(y(i,:)+0.5*h*k_2));
%     k_4 = fun((z(i)+h),(y(i,:)+k_3*h));
% 
%     y(i+1,:) = y(i,:) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  % main equation
% end
%%
