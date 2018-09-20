%% Correlations .

%% Obs 
for i =1:length(distance_class)
    for j =1:2
        for n = 1:500
            y = randsample(obs_pairs(j,i),obs_pairs(j,i),true);
            
            Obs_cor(n,j,i) = correlation_func(obs_sep(j,i).sep(y), ndays);
            
            Cobs(:,n,j,i) = Obs_cor(n,j,i).c; 
            robs(:,n,j,i) = Obs_cor(n,j,i).r; 
            
            % Mvel_obs(:,n,j,i) = Mvel_obs(:,n,j,i)/Mvel_obs(1,n,j,i);
        end
    end
end
%% model  
for i =1:length(distance_class)
    for j=1:2 
        Mod_cor(j,i) = correlation_func(mod_sep(j,i).sep, ndays) ;
        Cmod(:,j,i) = Mod_cor(j,i).c; 
        rmod(:,j,i) = Mod_cor(j,i).r; 
    end
end

%% 

[Cobsmean, Cobsci] = errbar4shaded(Cobs);

robsmean = squeeze(nanmean(robs,2));

%% 
% Figure at first depth 
close all 
clear h g
figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    
    h(i)= shadedErrorBar_semilogx(T+0.001, Cobsmean(:,1,i), Cobsci(:,:,1,i), ...
                    {'-','linewidth',3,'color',colors(i*2,:)},1)
                hold all
    g(i) = semilogx(T+0.001, Cmod(:,1,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end
legend([h(1).mainLine,h(2).mainLine, h(3).mainLine, g], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) 
axis([1 100 -0.2 1])

set(gca,'FontSize', 20)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$<\delta$ $\bf{V}$ $(t)$.$\delta$ $\bf{V_o}>/$ $|\delta \bf{V_o}|^2$', 'Interpreter','Latex')

%% separation on x axis 
% Figure at first depth 
close all 
clear h g
figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    
    h(i)= shadedErrorBar_semilogx(robsmean(:,1,i)/1e3, Cobsmean(:,1,i), Cobsci(:,:,1,i), ...
                    {'-','linewidth',3,'color',colors(i*2,:)},1)
                hold all
    g(i) = semilogx(rmod(:,1,i)/1e3, Cmod(:,1,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end
legend([h(1).mainLine,h(2).mainLine, h(3).mainLine, g], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) 
axis([10 300 -0.2 1])

set(gca,'FontSize', 20)
xlabel('$r$ (km)', 'Interpreter','Latex') 
ylabel('$<\delta$ $\bf{V}$ $(t)$.$\delta$ $\bf{V_o}>/$ $|\delta \bf{V_o}|^2$', 'Interpreter','Latex')

%% 
% Figure at second  depth 
close all 
clear h g
figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    
    h(i)= shadedErrorBar_semilogx(T+0.001, Cobsmean(:,2,i), Cobsci(:,:,2,i), ...
                    {'-','linewidth',3,'color',colors(i*2,:)},1)
                hold all
    g(i) = semilogx(T+0.001, Cmod(:,2,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end
legend([h(1).mainLine,h(2).mainLine, h(3).mainLine, g], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) 
axis([1 100 -0.2 1])

set(gca,'FontSize', 20)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('Vel Corr')

ylabel('$<$ $\bf{u}_1$ .$\bf{u}_1$ $>/(<$ $|\bf{u}_1|^2>$ $<|\bf{u}_2|^2>)^{1/2}$ ', 'Interpreter', 'Latex')
%ylabel('$<\delta$ $\bf{V}$ $(t)$.$\delta$ $\bf{V_o}>/$ $|\delta \bf{V_o}|^2$', 'Interpreter','Latex')

%% separation on x axis 
% Figure at first depth 
close all 
clear h g
figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    
    h(i)= shadedErrorBar_semilogx(robsmean(:,2,i)/1e3, Cobsmean(:,2,i), Cobsci(:,:,2,i), ...
                    {'-','linewidth',3,'color',colors(i*2,:)},1)
                hold all
    g(i) = semilogx(rmod(:,2,i)/1e3, Cmod(:,2,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end
legend([h(1).mainLine,h(2).mainLine, h(3).mainLine, g], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) 
axis([10 300 -0.2 1])

set(gca,'FontSize', 20)
xlabel('$r$ (km)', 'Interpreter','Latex') 
ylabel('$<\delta$ $\bf{V}$ $(t)$.$\delta$ $\bf{V_o}>/$ $|\delta \bf{V_o}|^2$', 'Interpreter','Latex')

%% Time Scale 

clear Tseries Ts_int Tscale_obs
%%
Tscale_obs = zeros(500, 2,6); 

for i =1:length(distance_class)
    for j =1:2
        for n = 1:500
            Tseries = Cobs(:,n,j,i);

            tax = linspace(0,99,1000); 

            Ts_int = interp1(T, Tseries, tax); 

            Tscale_obs(n,j,i) = tax(find(Ts_int<0.5,1));
        end
    end
end

for i =1:length(distance_class)
    for j =1:2

            Tseries = Cmod(:,j,i);

            tax = linspace(0,99,1000); 

            Ts_int = interp1(T, Tseries, tax); 

            Tscale_mod(j,i) = tax(find(Ts_int<0.5,1));
    end
end



%%
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

axis([0 70 0 70])
A = legend('Model 750m', 'Model 1500m', 'Obs 500-1000m', 'Obs 1000-1800m')
set(A, 'location', 'northwest', 'fontsize',16)
set(gca, 'Fontsize',20)
xlabel('D_o (km)')
ylabel('Mem. Time Scale (Days)')
saveas(gcf,'../figures/corr_time.eps', 'epsc')


%%
%%
Rscale_obs = zeros(500, 2,6); 

for i =1:length(distance_class)
    for j =1:2
        for n = 1:500
            Tseries = Cobs(:,n,j,i);
            rseries = robs(:,n,j,i)/1e3; 
            
            rax = linspace(10,150,1000); 

            Rs_int = interp1(rseries, Tseries, rax); 

            Rscale_obs(n,j,i) = rax(find(Rs_int<0.5,1));
        end
    end
end

%%
Rscale_obs_mean = squeeze(nanmean(Rscale_obs, 1)); 
for i =1:length(distance_class)
    for j =1:2
            Rscale_obs_ci(:,j,i) = prctile(Rscale_obs(:,j,i), [95, 5]);
    end
end
%%
Rscale_mod = zeros( 2,6); 

for i =1:length(distance_class)
    for j =1:2
        
            Tseries = Cmod(:,j,i);
            rseries = rmod(:,j,i)/1e3; 
            
            rax = linspace(10,150,1000); 

            Rs_int = interp1(rseries, Tseries, rax); 

            Rscale_mod(j,i) = rax(find(Rs_int<0.5,1));
        
    end
end

%%
%close all
dist_axis = mod_dist_ini(1,:)/1000;
figure('rend','painters','pos',[10 10 900 300])
plot(dist_axis, Rscale_mod,'*-','MarkerSize',10, 'Linewidth',3)
hold all 
errorbar(dist_axis, Rscale_obs_mean(1,:), ...
    squeeze(Rscale_obs_ci(2,1,:)), ...
    squeeze(Rscale_obs_ci(1,1,:)) ...
    ,'o-', 'MarkerSize',10, 'LineWidth',3)
errorbar(dist_axis, Rscale_obs_mean(2,:), ...
    squeeze(Rscale_obs_ci(2,2,:)), ...
    squeeze(Rscale_obs_ci(1,2,:)) ...
    ,'o-', 'MarkerSize',10, 'LineWidth',3)

%errorbar(dist_axis, Tscale_obs_mean(1,:), ...
%    squeeze(Tscale_obs_ci(2,1,:))'- Tscale_obs_mean(1,:), ...
%    squeeze(Tscale_obs_ci(1,1,:))'- Tscale_obs_mean(1,:) ...
%    ,'o-', 'MarkerSize',10, 'LineWidth',3)


%errorbar(dist_axis, Tscale_obs_mean(2,:), ...
%    squeeze(Tscale_obs_ci(2,2,:))'- Tscale_obs_mean(2,:), ...
%    squeeze(Tscale_obs_ci(1,2,:))'- Tscale_obs_mean(2,:) ...
%    ,'o-', 'MarkerSize',10, 'LineWidth',3)

%axis([0 70 0 70])
A = legend('Model 750m', 'Model 1500m', 'Obs 500-1000m', 'Obs 1000-1800m')
set(A, 'location', 'northwest', 'fontsize',16)
set(gca, 'Fontsize',20)
xlabel('D_o (km)')
ylabel('Mem. Time Scale (Days)')
%saveas(gcf,'../figures/corr_time.eps', 'epsc')

