function demi_vie = Demi_vie (X,Transfo)
% https://www.nndc.bnl.gov/nudat2/chartNuc.jsp
        Y = 365*24*60*60; % in [s]
        d = 24*60*60;     % in [s]
        h = 60*60;        % in [s]
        m = 60;           % in [s]
    switch X
        case 'U235'
            if (Transfo == 'Alpha')
                demi_vie=7.04 *10^8*Y;
            end
        case 'U238'
            if (Transfo == 'Alpha')
                demi_vie=4.468*10^9*Y;
            end
        case 'U239'
            if (Transfo == 'BetaMinus')
                demi_vie=23.45*m;
            end
        case 'Np239'
            if (Transfo == 'BetaMinus')
                demi_vie=2.356*d;
            end
        case 'Pu239'
            if (Transfo == 'Alpha')
                demi_vie=24110*Y;
            end
        case 'Xe135'
            if (Transfo == 'BetaMinus')
                demi_vie=9.14*h;
            end
        otherwise 
            disp('Element not in the list');    
      end
end
