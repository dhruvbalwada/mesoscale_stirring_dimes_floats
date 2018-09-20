%% Memory indexs

% obs
for i =1:length(distance_class)
    for j =1:2
        for n = 1:500
            y = randsample(obs_pairs(j,i),obs_pairs(j,i),true);
            [M_obs(:,n,j,i), Mvel_obs(:,n,j,i)] = memory_index(obs_sep(j,i).sep(y), ndays);
            %Mvel_obs(:,n,j,i) = Mvel_obs(:,n,j,i)/Mvel_obs(1,n,j,i);
        end
    end
end

%% mod
for i =1:length(distance_class)
    for j =1:2
         [M_mod(:,j,i), Mvel_mod(:,j,i)] = memory_index(mod_sep(j,i).sep, ndays);
          % Mvel_mod(:,j,i) = Mvel_mod(:,j,i)/Mvel_mod(1,j,i);          
    end
disp(i)
end

%% calc the means and std

Mobs_mean = squeeze(nanmean(M_obs,2));
Mobs_std = squeeze(nanstd(M_obs,0,2));

Mvelobs_mean = squeeze(nanmean(Mvel_obs,2));
Mvelobs_std = squeeze(nanstd(Mvel_obs,0,2));
%% percentiles as we don't expect the data to be gaussian
for i =1:length(distance_class)
    for j =1:2
        for k =1:ndays
            Mvelobsci(k,:,j,i) = prctile(Mvel_obs(k,:,j,i), [95, 5]);
        end
    end
end

% error bars for plotting 

for i =1:length(distance_class)
    for j =1:2
        Mvelobscir(:,1,j,i) = Mvelobsci(:,1,j,i) - Mvelobs_mean(:,j,i);
        Mvelobscir(:,2,j,i) = -Mvelobsci(:,2,j,i) + Mvelobs_mean(:,j,i);
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
    
    h(i)= shadedErrorBar_semilogx(T+0.001, Mvelobs_mean(:,1,i), Mvelobscir(:,:,1,i), ...
                    {'-','linewidth',3,'color',colors(i*2,:)},1)
                hold all
    g(i) = semilogx(T+0.001, Mvel_mod(:,1,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end
legend([h(1).mainLine,h(2).mainLine, h(3).mainLine, g], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) 
axis([1 100 -0.5 1])

set(gca,'FontSize', 20)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$<\delta$ $\bf{V}$ $(t)$.$\delta$ $\bf{V_o}>/$ $|\delta \bf{V_o}|^2$', 'Interpreter','Latex')

saveas(gcf,'../figures/memory_shallow.eps', 'epsc')
%%
% Figure at first depth 

figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    h(i)= shadedErrorBar_semilogx(T+0.001, Mvelobs_mean(:,2,i), Mvelobscir(:,:,2,i), ...
                    {'-','linewidth',3,'color',colors(i*2,:)},1)
                hold all
    g(i) = semilogx(T+0.001, Mvel_mod(:,2,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    i=i+1
    disp(k)
end
legend([h(1).mainLine,h(2).mainLine, h(3).mainLine, g], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) 
axis([1 100 -0.5 1])

set(gca,'FontSize', 20)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$<\delta$ $\bf{V}$ $(t)$.$\delta$ $\bf{V_o}>/$ $|\delta \bf{V_o}|^2$', 'Interpreter','Latex')

saveas(gcf,'../figures/memory_deep.eps', 'epsc')
%% Derive the 0.5 crossing time scale 

for i =1:length(distance_class)
    for j =1:2
        for n = 1:500
            Tseries = Mvel_obs(1:20,n,j,i);

            tax = linspace(0,19,10000); 

            Ts_int = interp1(T(1:20), Tseries, tax); 

            Tscale_obs(n,j,i) = tax(find(Ts_int<0.5,1));
        end
    end
end

for i =1:length(distance_class)
    for j =1:2

            Tseries = Mvel_mod(1:20,j,i);

            tax = linspace(0,19,1000); 

            Ts_int = interp1(T(1:20), Tseries, tax); 

            Tscale_mod(j,i) = tax(find(Ts_int<0.5,1));
    end
end

%% means and errors on time scales 

Tscale_obs_mean = squeeze(nanmean(Tscale_obs, 1)); 
for i =1:length(distance_class)
    for j =1:2
            Tscale_obs_ci(:,j,i) = prctile(Tscale_obs(:,j,i), [95, 5]);
    end
end

%% 
%close all
dist_axis = mod_dist_ini(1,:)/1000;
figure('rend','painters','pos',[10 10 900 300])
plot(dist_axis, Tscale_mod,'*-','MarkerSize',10, 'Linewidth',3)
hold all 
errorbar(dist_axis, Tscale_obs_mean(1,:), ...
    squeeze(Tscale_obs_ci(2,1,:)), ...
    squeeze(Tscale_obs_ci(1,1,:)) ...
    ,'o-', 'MarkerSize',10, 'LineWidth',3)
errorbar(dist_axis, Tscale_obs_mean(2,:), ...
    squeeze(Tscale_obs_ci(2,2,:)), ...
    squeeze(Tscale_obs_ci(1,2,:)) ...
    ,'o-', 'MarkerSize',10, 'LineWidth',3)

%errorbar(dist_axis, Tscale_obs_mean(1,:), ...
%    squeeze(Tscale_obs_ci(2,1,:))'- Tscale_obs_mean(1,:), ...
%    squeeze(Tscale_obs_ci(1,1,:))'- Tscale_obs_mean(1,:) ...
%    ,'o-', 'MarkerSize',10, 'LineWidth',3)


%errorbar(dist_axis, Tscale_obs_mean(2,:), ...
%    squeeze(Tscale_obs_ci(2,2,:))'- Tscale_obs_mean(2,:), ...
%    squeeze(Tscale_obs_ci(1,2,:))'- Tscale_obs_mean(2,:) ...
%    ,'o-', 'MarkerSize',10, 'LineWidth',3)

axis([0 70 0 18])
A = legend('Model 750m', 'Model 1500m', 'Obs 500-1000m', 'Obs 1000-1800m')
set(A, 'location', 'northwest', 'fontsize',16)
set(gca, 'Fontsize',20)
xlabel('D_o (km)')
ylabel('Mem. Time Scale (Days)')
saveas(gcf,'../figures/memory_time.eps', 'epsc')

%% 

