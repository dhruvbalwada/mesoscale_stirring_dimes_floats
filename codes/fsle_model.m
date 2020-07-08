% Dhruv Balwada
% 24 Septmber 2018

% sep is the collection of all possible pairs (no initial distance crit)
% that can be formed using the data.
% dist_bin - the distance bins
% limited to
% flag_1time - ???

function [fsle] = fsle_model(sep, dist_bin, flag_interp, flag_1time)

dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));

for i =1:length(dist_axis)
    
    lambda(i).ts = [];
    
    % loop over different pairs
    for j = 1:length(sep)
        
        
        % find id of the pairs in a particular geographical regime (lon
        % and pressure and within a certain distance from each other)
        id1 = find(sep(j).dist(1:end-1)<=dist_bin(i) & sep(j).dist(2:end)>=dist_bin(i));        

        id2 = find(sep(j).dist(1:end-1)<=dist_bin(i+1) & sep(j).dist(2:end)>=dist_bin(i+1))+1;      
        
        
        if ~isempty(id1) && ~isempty(id2)
            
            id22 = nan(length(id1),1);
            
            % for every id1 find the corresponding end of separation (if it
            % exists).
            for k =1:length(id1)
                id = find(id2>id1(k),1);
                if (~isempty(id))
                    id22(k) = id2(id);
                else
                    id22(k) = NaN;
                end
            end
            
            clear pair 
            m=1;
            
            
            for k = 1:length(id22)
                
                if isnan(id22(k))
                    continue
                end
                
                if flag_interp==0
                    %t1 = id1(k);
                    %t2 = id22(k);
                    t1 = sep(j).T(id1(k));
                    t2 = sep(j).T(id22(k));
                else
                    t1 = interp1([sep(j).dist(id1(k)) sep(j).dist(id1(k)+1)], [sep(j).T(id1(k)) sep(j).T(id1(k)+1)], dist_bin(i));
                    t2 = interp1([sep(j).dist(id22(k)-1) sep(j).dist(id22(k))], [sep(j).T(id22(k)-1)  sep(j).T(id22(k))], dist_bin(i+1));
                end
                
                if (t2>t1)
                    pair(m) = t2 - t1;
                end
                
                % if this flag is set to 1, then the fastest time to make the separation is only considered.
                if flag_1time==1
                    continue
                end
                
                if k<length(id22)
                    while id22(k)==id22(k+1)
                        k=k+1;
                        t1 = interp1([sep(j).dist(id1(k)) sep(j).dist(id1(k)+1)], [sep(j).T(id1(k)) sep(j).T(id1(k)+1)], dist_bin(i));
                        t2 = interp1([sep(j).dist(id22(k)-1) sep(j).dist(id22(k))], [sep(j).T(id22(k)-1)  sep(j).T(id22(k))], dist_bin(i+1));
                        
                        if t2>t1
                            pair(m) = t2 - t1;
                        end
                        
                        if k==length(id22)
                            break
                        end
                    end
                end
                m = m+1;
            end
            
            
            if exist('pair')
                lambda(i).ts = [lambda(i).ts log(dist_bin(i+1)/dist_bin(i))./pair];
            end
            
        end
    end
    %disp(i)
    
end

for i = 1:length(dist_axis)
    fsle(i) = nanmean(lambda(i).ts);
end


end