% Dhruv Balwada
% 24 September 2018

function [S2llmean, S2llci] = S2ll(sep, dist_bin, diff_pres, plevel)

dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));

for i =1:length(dist_axis)
    dull = [];
    % loop over different pairs
    for j = 1:length(sep)
        dull_temp = [];
        % find id of the pairs in a particular geographical regime (lon
        % and pressure and within a certain distance from each other)
        
        id = find(sep(j).dist<dist_bin(i+1) & sep(j).dist>=dist_bin(i) ...
            & sep(j).dP<diff_pres & ...
            sep(j).P1>plevel(1) & sep(j).P1<plevel(2));

        
        ids2keep = find(id~=1117);
        %  id = id(ids2keep);
        
        % loop over the different pairs that lie in the range
        for k =1:length(id)
            % components
            rx(k) = (sep(j).X1(id(k)) - sep(j).X2(id(k)))*cosd(0.5*(sep(j).Y1(id(k))+sep(j).Y2(id(k))));
            ry(k) = (sep(j).Y1(id(k)) - sep(j).Y2(id(k)));
            magr(k) = sqrt(rx(k).^2+ry(k).^2);
            % normalize to unit vectors
            rx(k) = rx(k)/magr(k); ry(k) = ry(k)/magr(k);
            
            % components of velocity differences
            dux(k) = (sep(j).U1(id(k))-sep(j).U2(id(k)));
            duy(k) = (sep(j).V1(id(k))-sep(j).V2(id(k)));
            
            % convert to longitudnal and
            dull_temp(k) = dux(k)*rx(k) + duy(k)*ry(k);
            
        end
        if ~isempty(dull_temp)
            dull = [dull; dull_temp'];
            
        end
    end
    
    struct_pairs(i).dull = dull;
    
    disp(i)
end

for i = 1:length(dist_axis)
    for n=1:500 
        len = length(struct_pairs(i).dull); 
        y = randsample(len, len, true); 
        
        if len>100
            S2llbs(i,n) = nanmean(struct_pairs(i).dull(y).^2);
        else 
            S2llbs(i,n) = nan;
        end
    end
end

S2llmean = nanmean(S2llbs, 2); 

for i = 1:length(dist_axis) 
    S2lltemp(i, :) = prctile(S2llbs(i,:), [95, 5]); 
end

S2llci(:,1) = 0.434*(S2lltemp(:,1) - S2llmean)./S2llmean; 
S2llci(:,2) = 0.434*(-S2lltemp(:,2) + S2llmean)./S2llmean;

end