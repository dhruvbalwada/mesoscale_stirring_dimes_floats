% Modified on 27 January 2020
% Dhruv Balwada
% 
% rel_diff.m
% 
% Calculate the relative diffusivity in the simple way 
% 
% Krel = 1/2 * d<Y^2>/dt

function [diff_struct] = rel_diff(dispersion, ndays, smdays) 
    
diff = nan*ones(ndays,1); 
    r = nan*ones(ndays,1); 
    
    diff(1:end-smdays) = (dispersion(1+smdays:end) - dispersion(1:end-smdays))/24/3600/2/smdays; 
    
    r = dispersion.^0.5; 
    
   diff_struct.diff =  diff; 
   diff_struct.r =  r; 
end