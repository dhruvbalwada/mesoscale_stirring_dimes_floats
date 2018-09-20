% Modified on 30 August 2018
% Dhruv Balwada
% 
% rel_diff.m
% 
% Calculate the relative diffusivity in the simple way 
% 
% Krel = 1/2 * d<Y^2>/dt

function [diff] = rel_diff(dispersion, dt, ndays) 
    diff = nan*ones(ndays,1); 
    
    dtm1 = dt - 1;    
    
    diff(1+dtm1:end-dtm1) = (dispersion(1+dt:end) - dispersion(1:end-dt))/24/3600/2/dt; % center difference
    
end