% function to do the PDF fitting to the displacement PDFS
% do a least squares fit here vs the point fit that was done before.

function [fitting] = dispersion_fitting_log(dispersion,T,a)

%global Ro
% r_o. Initial separation and time when the separation has increased to
% a times that (in a least squares sense).

Ro = sqrt(dispersion(1)); 
% id when the
id = find(dispersion>=a*(Ro^2),1);
tstar = T(id);
%% Matching the Lungdren Distribution (least squares)

Tlguess = (8*tstar)/(log(a));

% function for the dispersion of non-local dispersion
flung = @(Tlung,xdata)log((Ro^2)*exp(8*xdata/Tlung));

Tlung = lsqcurvefit(flung, Tlguess, T(1:id), log(dispersion(1:id)));

%disp(flung(Tlung, T(1)))
%% Matching to the Richardson Distribution (least squares)

bguess = 27*Ro^(2/3)/(20*tstar);

% function for the dispersion of the local dispersion
frich = @(bguess,xdata) log( factorial(5)/2* ...
    ((4*bguess*xdata/9).^3).* ...
    exp(-9*Ro^(2/3)/4/bguess./xdata).* ...
    hypergeom(6,3,9*Ro^(2/3)/4/bguess./xdata));

beta = lsqcurvefit(frich,bguess, T(2:id), log(dispersion(2:id)));
%disp(frich(beta, T(2))) 
%% Matching the Rayleigh Distribution
%

% function of the dispersion for the diffusive dispersion

Kguess = dispersion(id)/4/tstar;

fdiff = @(kappa, xdata) log(Ro^2 + 4*kappa*xdata);

kappa = lsqcurvefit(fdiff, Kguess, T(1:id), log(dispersion(1:id)));

%disp(fdiff(kappa,T(1)))

fitting.tstar = tstar;
fitting.Tlung = Tlung;
fitting.beta = beta;
fitting.kappa = kappa;
fitting.Ro = Ro;

end