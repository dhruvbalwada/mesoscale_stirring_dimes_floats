% function to do the PDF fitting to the displacement PDFS
% do a least squares fit here vs the point fit that was done before.

function [fitting] = dispersion_fitting_overtime(dispersion,T,Tmin, Tmax)

%global Ro
% r_o. Initial separation and time when the separation has increased to
% a times that (in a least squares sense).

Ro = sqrt(dispersion(1));

Tfit = Tmin:Tmax;

dispT0 = Ro^2;

for i = 1:length(Tfit)
    
    id = find(T==Tfit(i));
    
    tstar = Tfit(i);
    dispT = dispersion(id);
    a = dispT/dispT0;

    aT(i) = a; 
    
    %% Matching the Lungdren Distribution (least squares)
    
    Tlguess = (8*tstar)/(log(a));
    
    % function for the dispersion of non-local dispersion
    flung = @(Tlung,xdata)(Ro^2)*exp(8*xdata/Tlung);
    
    Tlung(i) = lsqcurvefit(flung, Tlguess, T(1:id), dispersion(1:id));
    
    %disp(flung(Tlung, T(1)))
    %% Matching to the Richardson Distribution (least squares)
    
    bguess = 27*Ro^(2/3)/(20*tstar);
    
    % function for the dispersion of the local dispersion
    frich = @(bguess,xdata) factorial(5)/2* ...
        ((4*bguess*xdata/9).^3).* ...
        exp(-9*Ro^(2/3)/4/bguess./xdata).* ...
        hypergeom(6,3,9*Ro^(2/3)/4/bguess./xdata);
    
    beta(i) = lsqcurvefit(frich,bguess, T(2:id), dispersion(2:id));
    %disp(frich(beta, T(2)))
    %% Matching the Rayleigh Distribution
    %
    
    % function of the dispersion for the diffusive dispersion
    
    %Kguess = dispersion(id)/4/tstar;
    
    %fdiff = @(kappa, xdata) Ro^2 + 4*kappa*xdata;
    
    %kappa(i) = lsqcurvefit(fdiff, Kguess, T(1:id), dispersion(1:id));
    
    %disp(fdiff(kappa,T(1)))
end

fitting.Tlung = Tlung;
fitting.beta = beta;
fitting.Ro = Ro;
fitting.aT = aT; 

end

