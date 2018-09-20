% Dhruv Balwada
% 18 September 2018

for i =1:length(distance_class)
    for j =1:2
        for n = 1:500
            y = randsample(obs_pairs(j,i),obs_pairs(j,i),true);
            C_obs(:,n,j,i) = chance_correlation(obs_sep(j,i).sep(y), ndays);
            % Mvel_obs(:,n,j,i) = Mvel_obs(:,n,j,i)/Mvel_obs(1,n,j,i);
        end
    end
    endl

%% 
for i =1:length(distance_class)
    for j =1:2
        
        
            C_mod(:,j,i) = chance_correlation(mod_sep(j,i).sep, ndays);
            % Mvel_obs(:,n,j,i) = Mvel_obs(:,n,j,i)/Mvel_obs(1,n,j,i);
        
    end
end

%%

Cobs_mean = squeeze(nanmean(C_obs,2));

%% percentiles as we don't expect the data to be gaussian
for i =1:length(distance_class)
    for j =1:2
        for k =1:ndays
            Cobsci(k,:,j,i) = prctile(C_obs(k,:,j,i), [95, 5]);
        end
    end
end

% error bars for plotting 

for i =1:length(distance_class)
    for j =1:2
        Cobscir(:,1,j,i) = Cobsci(:,1,j,i) - Cobs_mean(:,j,i);
        Cobscir(:,2,j,i) = -Cobsci(:,2,j,i) + Cobs_mean(:,j,i);
    end
end

%%
colors =get(gca,'colororder');
%% 
% Figure at first depth 
close all 
clear h g
figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    
    h(i)= shadedErrorBar_semilogx(T+0.001, Cobs_mean(:,2,i), Cobscir(:,:,2,i), ...
                    {'-','linewidth',3,'color',colors(i*2,:)},1)
                hold all
    g(i) = semilogx(T+0.001, C_mod(:,2,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end
legend([h(1).mainLine,h(2).mainLine, h(3).mainLine, g], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) 
axis([1 100 -0.5 0.5])

set(gca,'FontSize', 20)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$<\delta$ $\bf{V}$ $(t)$.$\delta$ $\bf{V_o}>/$ $|\delta \bf{V_o}|^2$', 'Interpreter','Latex')

%saveas(gcf,'../figures/chance_corr.eps', 'epsc')
