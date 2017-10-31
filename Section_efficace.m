function [sigma] = Section_efficace(X,Transfo,n_eV,User_Adress)
% X, espèce chimique
% Transfo, Fission/Capture
% n_eV, energie du neutron incident
% User_Adress, adresse de la DB
    item = 0;
    X = 'Pu239'
    Transfo = 'Fission'
    switch X
        case 'Pu239'
            item = item + 1;
        otherwise
            item = item + 0;
        
    end
    switch Transfo
        case 'Capture'
            item = item +10;
        case 'Fission'
            item = item +20
        otherwise
            item = item + 0;
    end

    if item == 21
        A = load('data/Pu239F.txt');
        loglog(A(:,1),A(:,2));
        

    end    
end
