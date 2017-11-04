function [sigma] = Section_efficace(X,Transfo,n_eV,User_Adress)
% X, espèce chimique
% Transfo, Fission/Capture
% n_eV, energie du neutron incident
% User_Adress, adresse de la DB
% 

    filename = sprintf('data/%s_%s.txt',X,Transfo);
    DATA = load(filename);    
    %% Chercheur et Interpolateur
    flag = 0;
    i = 1;
    while (i < length(DATA) && flag == 0) 
        if  (DATA(i,1)<= n_eV  && n_eV <= DATA(i+1,1))
            a = (DATA(i+1,2)-DATA(i,2))/(DATA(i+1,1)-DATA(i,1));
            b = (DATA(i+1,2)*DATA(i,1) - DATA(i,2)*DATA(i+1,1))/(DATA(i,1)-DATA(i+1,1));
            sigma = a*n_eV + b
            flag = 1;
        end
        i = i+1;
    end


 
end

%%Voilà ce que j'ai fais, on peut modifier pour intégrer l'idée de Son :) 

function [sigma] = Section_efficace(X,Transfo,n_eV,User_Adress)

User_Adress = strcat(User_Adress,X);

if Transfo == 'Capture'
    if X == 'Pu9'
        disp('Capture Pu9 pas possible');
        return;
    else
        User_Adress = strcat(User_Adress,'_Capture.txt');
    end
else
    User_Adress = strcat(User_Adress,'_Fission.txt');
end

tableau = load(User_Adress);
energy = abs(n_eV - tableau(:,1));
section = tableau(:,2);

[minEnergy index] = min(energy(:,1));
sigma = section(index);

end
