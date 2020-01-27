% Dhruv Balwada
% 14 Jan 2020

function [npairs_dist] = npairs_distbins(sep, dist_bin, diff_pres, plevel)

dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));

npairs_dist = 0*dist_axis; 

for i =1:length(dist_axis)
    % loop over different pairs
    
    for j = 1:length(sep)
      
        % find id of the pairs in a particular geographical regime (lon
        % and pressure and within a certain distance from each other)
        
        id = find(sep(j).dist<dist_bin(i+1) & sep(j).dist>=dist_bin(i) ...
            & sep(j).dP<diff_pres & ...
            sep(j).P1>plevel(1) & sep(j).P1<plevel(2));

        npairs_dist(i) = npairs_dist(i) + length(id); 
    
    end
    
end

end