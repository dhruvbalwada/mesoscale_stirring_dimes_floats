% Dhruv Balwada
% 3 September 2018 

% Calculate the temporal structure function 
% <v^2(t)> = <dr/dt^2>=<((u_i - u_j).r)^2>
% Averaging done over all pairs that are initially separated by specified
% separation distance 

function [S2_1, S2_2] = time_structure_function(sep, ndays)
    
    % method 1
    for i = 1:length(sep)
        sep(i).drdt = diff(sep(i).dist)/24/3600;
        
        drdt(1:ndays,i) = sep(i).drdt(1:ndays); 
    end
    
    S2_1 = nanmean(drdt.^2, 2);
        
    % method 2
    for i = 1:length(sep)
        % dull_temp(k) = dux(k)*rx(k) + duy(k)*ry(k);

        sep(i).dull = (sep(i).dU.*sep(i).X + sep(i).dV.*sep(i).Y)./sep(i).dist; 
        
        dull(1:ndays,i) = sep(i).dull(1:ndays); 
    end
    
    S2_2 = nanmean(dull.^2, 2);
end