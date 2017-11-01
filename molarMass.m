function M = molarMass(X)


% Input : [string] : nuclear element :
% 'U235','U238','U239','Np239','Pu239','Xe135'.

%Output : [double] : Molar mass of the element.
%v1.0 31/10/17
%%

    switch X
        case 'U235'
            M=235.043931368;
        case 'U238'
            M=238.050789466;
        case 'U239'
            M=239.054294518;
        case 'Np239'
            M=239.052940487;
        case 'Pu239'
            M=239.052164844;
        case 'Xe135'
            M=134.907226844;
        otherwise 
            disp('Element not in the list');    
    end
end