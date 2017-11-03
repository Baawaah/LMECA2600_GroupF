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
