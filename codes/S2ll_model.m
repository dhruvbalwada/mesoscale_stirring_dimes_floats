% Dhruv Balwada
% 24 September 2018

function [S2ll] = S2ll_model(sep, dist_bin)

dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));

for i =1:length(dist_axis)
    dull = [];
    % loop over different pairs
    for j = 1:length(sep)
        dull_temp = [];
        % find id of the pairs in a particular geographical regime (lon
        % and pressure and within a certain distance from each other)
        
        id = find(sep(j).dist<dist_bin(i+1) & sep(j).dist>=dist_bin(i));
        
        %         id = find(sep(j).dist<dist_bin(i+1) & sep(j).dist>=dist_bin(i) ...
        %             & abs(sep(j).T1-sep(j).T2)<diff_temp & sep(j).X1<-70 & sep(j).X2<-70 & ...
        %             sep(j).P1>plevel(1) & sep(j).P1<plevel(2));
        
        %ids2keep = find(id~=1117);
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
    
    %disp(i)
end

for i = 1:length(dist_axis)
   
            S2ll(i) = nanmean(struct_pairs(i).dull.^2);
end


end


