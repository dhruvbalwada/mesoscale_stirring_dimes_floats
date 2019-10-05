% Dhruv Balwada
% 24 September 2018

function [S3_vec] = S3(sep, dist_bin, diff_pres, plevel)

dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));

for i =1:length(dist_axis)
    dull = [];
    dutt = [];
    % loop over different pairs
    for j = 1:length(sep)
        dull_temp = [];
        dutt_temp = []; 
        % find id of the pairs in a particular geographical regime (lon
        % and pressure and within a certain distance from each other)
        
        id = find(sep(j).dist<dist_bin(i+1) & sep(j).dist>=dist_bin(i) ...
            & sep(j).dP<diff_pres & ...
            sep(j).P1>plevel(1) & sep(j).P1<plevel(2));
        
        %         id = find(sep(j).dist<dist_bin(i+1) & sep(j).dist>=dist_bin(i) ...
        %             & abs(sep(j).T1-sep(j).T2)<diff_temp & sep(j).X1<-70 & sep(j).X2<-70 & ...
        %             sep(j).P1>plevel(1) & sep(j).P1<plevel(2));
        
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
            
            % convert to longitudnal 
            dull_temp(k) = dux(k)*rx(k) + duy(k)*ry(k);
            dutt_temp(k) = duy(k)*rx(k) - dux(k)*ry(k);
            
        end
        if ~isempty(dull_temp)
            dull = [dull; dull_temp'];
            dutt = [dutt; dutt_temp'];
        end
    end
    
    struct_pairs(i).dull = dull;
    struct_pairs(i).dutt = dutt;
    
    disp(i)
end

for i = 1:length(dist_axis)
    S3_vec.S3LLL(i) = nanmean(struct_pairs(i).dull.^3); 
    S3_vec.S3LLT(i) = nanmean(struct_pairs(i).dull.^2.*struct_pairs(i).dutt); 
    S3_vec.S3LTT(i) = nanmean(struct_pairs(i).dutt.^2.*struct_pairs(i).dull); 
    S3_vec.S3TTT(i) = nanmean(struct_pairs(i).dutt.^3); 
end

% commented code is below for calculating error bars using bootstrapping. 
% for i = 1:length(dist_axis)
%     for n=1:500 
%         len = length(struct_pairs(i).dull); 
%         y = randsample(len, len, true); 
%         
%         if len>100
%             S3llbs(i,n) = nanmean(struct_pairs(i).dull(y).^3);
%         else 
%             S3llbs(i,n) = nan;
%         end
%     end
% end
% 
% S3llmean = nanmean(S3llbs, 2); 
% 
% for i = 1:length(dist_axis) 
%     S3lltemp(i, :) = prctile(S3llbs(i,:), [95, 5]); 
% end
% 
% S3llci(:,1) = 0.434*(S3lltemp(:,1) - S3llmean)./S3llmean; 
% S3llic(:,2) = 0.434*(-S3lltemp(:,2) + S3llmean)./S3llmean;

end


