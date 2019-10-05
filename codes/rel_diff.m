% Modified on 30 August 2018
% Dhruv Balwada
% 
% rel_diff.m
% 
% Calculate the relative diffusivity in the simple way 
% 
% Krel = 1/2 * d<Y^2>/dt

function [diff] = rel_diff(dispersion, ndays) 
    diff = nan*ones(ndays,1); 
    
    
    
    diff(2:end-1) = (dispersion(3:end) - dispersion(1:end-2))/24/3600/2/2; % center difference
    
end