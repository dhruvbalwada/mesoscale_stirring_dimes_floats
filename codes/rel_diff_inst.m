% Dhruv Balwada
% 24 September 2018

function [rdmean_struct] = rel_diff_inst(sep, dist_bin, diff_pres, plevel, smdays)

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
        
        id = find(sep(j).dist<dist_bin(i+1) & sep(j).dist>=dist_bin(i) ...
            & sep(j).dP<diff_pres & ...
            sep(j).P1>plevel(1) & sep(j).P1<plevel(2));

        
        ids2keep = find(id~=1117);
        %  id = id(ids2keep);

        % smoothing of fields before calc differences 
%         X = nanmoving_average(sep(j).X,2);
%         Y = nanmoving_average(sep(j).Y,2);
%         dist = nanmoving_average(sep(j).dist,2);

        
        % loop over the different pairs that lie in the range
        for k =1:length(id)
            
            rdX_temp(k) = (sep(j).X(id(k)+smdays).^2 - sep(j).X(id(k)).^2 )/24/3600/smdays/2;
            rdY_temp(k) = (sep(j).Y(id(k)+smdays).^2 - sep(j).Y(id(k)).^2 )/24/3600/smdays/2;
            rd_temp(k)  = (sep(j).dist(id(k)+smdays).^2 - sep(j).dist(id(k)).^2 )/24/3600/smdays/2;
            
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
    
    %disp(i)
end


% Error bar calculations
for i = 1:length(dist_axis)
    for n=1:500 
        len = length(rdpairs(i).rd); 
        y = randsample(len, len, true); 
        
        if len>100
            rd_bs(i,n) = nanmean(rdpairs(i).rd(y));
            rdX_bs(i,n) = nanmean(rdpairs(i).rdX(y));
            rdY_bs(i,n) = nanmean(rdpairs(i).rdY(y));
        else 
            rd_bs(i,n) = nan;
            rdX_bs(i,n) = nan;
            rdY_bs(i,n) = nan;
        end
    end
end

rdmean = nanmean(rd_bs, 2); 
rdXmean = nanmean(rdX_bs, 2); 
rdYmean = nanmean(rdY_bs, 2); 

rdmean_struct.rdmean = rdmean;
rdmean_struct.rdXmean = rdXmean;
rdmean_struct.rdYmean = rdYmean;


for i = 1:length(dist_axis) 
    rdtemp(i, :) = prctile(rd_bs(i,:), [95, 5]);
    rdXtemp(i, :) = prctile(rdX_bs(i,:), [95, 5]);
    rdYtemp(i, :) = prctile(rdY_bs(i,:), [95, 5]);
end

rdci(:,1) = 0.434*(rdtemp(:,1) - rdmean)./rdmean; 
rdci(:,2) = 0.434*(-rdtemp(:,2) + rdmean)./rdmean;

rdXci(:,1) = 0.434*(rdXtemp(:,1) - rdXmean)./rdXmean; 
rdXci(:,2) = 0.434*(-rdXtemp(:,2) + rdXmean)./rdXmean;

rdYci(:,1) = 0.434*(rdYtemp(:,1) - rdYmean)./rdYmean; 
rdYci(:,2) = 0.434*(-rdYtemp(:,2) + rdYmean)./rdYmean;

rdmean_struct.rdci = rdci;
rdmean_struct.rdXci = rdXci;
rdmean_struct.rdYci = rdYci;

end