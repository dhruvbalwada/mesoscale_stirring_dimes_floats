% Dhruv Balwada
% 18 September 2018 

% model 
for i =1:length(distance_class)
    for j =1:2        
            
            Mod_disp = rel_disp(mod_sep(j,i).sep, ndays, 1);

    end
end


%% obs
for i =1:length(distance_class)
    for j =1:2
        for n = 1:500
            y = randsample(obs_pairs(j,i),obs_pairs(j,i),true);
            Obsbs = rel_disp(obs_sep(j,i).sep(y), ndays, 1);
        end
    end
end

%% 

%%
Obs_dispnorm = squeeze(nanmean(Obs_dispnormbs,2));

%%
for i =1:length(distance_class)
    for j =1:2
        for n =1:ndays
            Obs_dispnormci(:,n,j,i) = prctile(Obs_dispnormbs(n,:,j,i), [95, 5]);   
        end
    end
end

%%
for i =1:length(distance_class)
    for j =1:2
        
            Obs_dispnormcir(1,:,j,i) = 0.434*(Obs_dispnormci(1,:,j,i)' - Obs_dispnorm(:,j,i))./ ...                
                                        Obs_dispnorm(:,j,i);   
            Obs_dispnormcir(2,:,j,i) = 0.434*(-Obs_dispnormci(1,:,j,i)' + Obs_dispnorm(:,j,i))./ ...                
                                        Obs_dispnorm(:,j,i);   
        
    end
end

%% 

colors = get(gca, 'ColorOrder')
%% 
close all 
figure 
for i =1:2:length(distance_class)
    for j =1
        semilogx(T, squeeze(Obs_dispnorm(:,j,i))./T'.^2,'o-', 'color',colors(i,:))
        hold all
        semilogx(T, squeeze(Mod_dispnorm(:,j,i))./T'.^2,'--',  'color',colors(i,:))
        
    end
end
%loglog(T, 10^7*T.^2, 'linewidth',3)


%% 
close all 
figure 
for i =1:2:length(distance_class)
    for j =1
        semilogx(T, squeeze(Obs_dispnorm(:,j,i))./T'.^3,'o-', 'color',colors(i,:))
        hold all
        semilogx(T, squeeze(Mod_dispnorm(:,j,i))./T'.^3,'--',  'color',colors(i,:))
    end
end


%% 
close all 
figure 
for i =1:2:length(distance_class)
    for j =1
        semilogy(T, squeeze(Obs_dispnorm(:,j,i)),'o-', 'color',colors(i,:))
        hold all
        semilogy(T, squeeze(Mod_dispnorm(:,j,i)),'--',  'color',colors(i,:))
    end
end


%% 
figure 
for i =1:2:length(distance_class)
    for j =2
        loglog(T, Mod_dispnorm(:,j,i),'o-')
hold all
    end
end
loglog(T, 10^7*T.^2, 'linewidth',3)