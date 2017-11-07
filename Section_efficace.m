function [sigma] = Section_efficace(X,Transfo,n_eV,User_Adress)
% X, espèce chimique
% Transfo, Fission/Capture
% n_eV, energie du neutron incident
% User_Adress, adresse de la DB
% 


if Transfo == 'Capture'
    if X == 'Pu9'
        disp('Capture Pu9 pas possible');
        return;
    end
end

    filename = sprintf('/%s_%s.txt',X,Transfo);
    filename = strcat(User_Adress,filename);
    DATA = load(filename);    
    %% Chercheur et Interpolateur
    flag = 0;
    i = 1; 
        while (i < length(DATA) && flag == 0) 
            if  (DATA(i,1)*10^6<= n_eV  && n_eV <= DATA(i+1,1)*10^6)
                a = (DATA(i+1,2)*10^6-DATA(i,2)*10^6)/(DATA(i+1,1)*10^6-DATA(i,1)*10^6);
                b = (DATA(i+1,2)*10^6*DATA(i,1)*10^6 - DATA(i,2)*10^6*DATA(i+1,1)*10^6)/(DATA(i,1)*10^6-DATA(i+1,1)*10^6);
                sigma = (a*n_eV + b)/10^6;
                flag = 1;
            end
            i = i+1;
        end
end


