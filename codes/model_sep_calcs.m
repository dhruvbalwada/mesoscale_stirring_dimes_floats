for dlvl = 1:2
    
    for d =1:length(distance_class)
        
        dist1 = distance_class(d).dist(1)/1e3;
        dist2 = distance_class(d).dist(2)/1e3;
        
        
        tj =[];
        tk =[];
        tdist = [];
        
        for i=1:12
            ids = (i-1)*100+1;
            ide = i*100;
            [dist, nj, nk, num] = pairs_identification(mod_traj(dlvl).X(:,ids:ide), mod_traj(dlvl).Y(:,ids:ide), dist1,dist2, 'y');
            tj = [tj nj+(i-1)*100];
            tk = [tk nk+(i-1)*100];
            tdist = [tdist dist];
        end

        %
        cor=(1.1120e5)^2;
        
        
        for i=1:length(tdist)
            x1 = mod_traj(dlvl).X(:,tj(i));
            x2 = mod_traj(dlvl).X(:,tk(i));
            y1 = mod_traj(dlvl).Y(:,tj(i));
            y2 = mod_traj(dlvl).Y(:,tk(i));
            
            separation=sqrt((cor)*(((x1-x2).*cosd((y1+y2)/2)).^2 + (y1-y2).^2));
            
            u1 = mod_traj(dlvl).U(:,tj(i));
            u2 = mod_traj(dlvl).U(:,tk(i));
            v1 = mod_traj(dlvl).V(:,tj(i));
            v2 = mod_traj(dlvl).V(:,tk(i));
            
            sep_final(i).X = (x1 - x2).*cosd((y1+y2)/2)*sqrt(cor) ;
            sep_final(i).Y = (y1 - y2)*sqrt(cor);
            sep_final(i).dist= separation;
            
            sep_final(i).dU = u1 - u2;
            sep_final(i).dV = v1 - v2;
            
            sep_final(i).X1 = x1;
            sep_final(i).Y1 = y1;
            sep_final(i).X2 = x2;
            sep_final(i).Y2 = y2;
            
            sep_final(i).U1 = u1;
            sep_final(i).V1 = v1;
            sep_final(i).U2 = u2;
            sep_final(i).V2 = v2;
        end
        
        
        mod_sep(dlvl,d).sep = sep_final;
        mod_sep(dlvl,d).depth = depth(dlvl);
        
        clear sep_final
    end
end