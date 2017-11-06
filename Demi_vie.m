function demi_vie = Demi_vie (X,Transfo)
% https://www.nndc.bnl.gov/nudat2/chartNuc.jsp
    switch X
        Y = 365*24*60*60; % in [s]
        d = 24*60*60;     % in [s]
        h = 60*60;        % in [s]
        m = 60;           % in [s]
        case 'U235'
            demi_vie=7.04 *10^8*Y;
        case 'U238'
            demi_vie=4.468*10^9*Y;
        case 'U239'
            demi_vie=     23.45*m;
        case 'Np239'
            demi_vie=     2.356*d;
        case 'Pu239'
            demi_vie=     24110*Y;
        case 'Xe135'
            demi_vie=      9.14*h;
        otherwise 
            disp('Element not in the list');    
    end
end