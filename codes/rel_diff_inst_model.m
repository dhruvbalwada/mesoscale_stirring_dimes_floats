% Dhruv Balwada
% 24 September 2018

function [rd_struct] = rel_diff_inst_model(sep, dist_bin,smdays)

dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));

for i =1:length(dist_axis)
    rd = [];
    rdX = [];
    rdY = [];
    
    % loop over different pairs
    for j = 1:length(sep)
        rdX_temp = [];
        rdY_temp = [];
        rd_temp = [];
        % find id of the pairs in a particular geographical regime (lon
        % and pressure and within a certain distance from each other)
        
        id = find(sep(j).dist<dist_bin(i+1) & sep(j).dist>=dist_bin(i));
        
        %         id = find(sep(j).dist<dist_bin(i+1) & sep(j).dist>=dist_bin(i) ...
        %             & abs(sep(j).T1-sep(j).T2)<diff_temp & sep(j).X1<-70 & sep(j).X2<-70 & ...
        %             sep(j).P1>plevel(1) & sep(j).P1<plevel(2));
        
        ids2keep = find(id<=734);
        id = id(ids2keep);
        
        % loop over the different pairs that lie in the range
        for k =1:length(id)
            % components
            rdX_temp(k) = (sep(j).X(id(k)+smdays).^2 - sep(j).X(id(k)).^2 )/24/3600/2/smdays;
            rdY_temp(k) = (sep(j).Y(id(k)+smdays).^2 - sep(j).Y(id(k)).^2 )/24/3600/2/smdays;
            rd_temp(k)  = (sep(j).dist(id(k)+smdays).^2 - sep(j).dist(id(k)).^2 )/24/3600/2/smdays;
            
        end
        if ~isempty(rd_temp)
            rdX = [rdX; rdX_temp'];
            rdY = [rdY; rdY_temp'];
            rd  = [rd; rd_temp'];
        end
    end
    
    rdpairs(i).rd = rd;
    rdpairs(i).rdX = rdX;
    rdpairs(i).rdY = rdY;
    
    disp(i)
end


for i = 1:length(dist_axis)
    rdmean(i) = nanmean(rdpairs(i).rd);
    rdXmean(i) = nanmean(rdpairs(i).rdX);
    rdYmean(i) = nanmean(rdpairs(i).rdY);
           %S2ll(i) = nanmean(struct_pairs(i).dull.^2);
end

rd_struct.rdmean = rdmean;
rd_struct.rdXmean = rdXmean;
rd_struct.rdYmean = rdYmean;

end


