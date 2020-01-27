function [ moment ] = moments_separation( seps, n )
% Estimate the nth moment of the separation curves: <r^n>

% we do the sqrt because the disps are separations^2.
moment = nanmean( seps'.^n, 1); 


end

