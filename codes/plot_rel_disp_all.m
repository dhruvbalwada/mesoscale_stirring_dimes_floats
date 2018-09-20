% Dhruv Balwada
% 18 September 2018 

% model 
for i =1:length(distance_class)
    for j =1:2        
            
            Mod_disp(j,i) = rel_disp(mod_sep(j,i).sep, ndays, 1);
            
    end
end

%% obs
for i =1:length(distance_class)
    for j =1:2
        for n = 1:500
            y = randsample(obs_pairs(j,i),obs_pairs(j,i),true);
            Obsbs(n,j,i) = rel_disp(obs_sep(j,i).sep(y), ndays, 1);
        end
    end
end

%% Isotropy 
plot_isotropy_all

%% Relative dispersion 

for i =1:length(distance_class)
    for j =1:2   
        Moddisp(:,j,i) = Mod_disp(j,i).avdisp; 
        Moddisp_zon(:,j,i) = Mod_disp(j,i).avdispzon;
        Moddisp_mer(:,j,i) = Mod_disp(j,i).avdispmer;
    end
end

for i =1:length(distance_class)
    for j =1:2   
        for n=1:500 
            Obsdisp(:,n,j,i) = Obsbs(n,j,i).avdisp;
            Obsdisp_zon(:,n,j,i) = Obsbs(n,j,i).avdispzon;
            Obsdisp_mer(:,n,j,i) = Obsbs(n,j,i).avdispmer;
        end
    end
end

%% 
[Obsdisp_mean, Obsdisp_ci] = errbar4shaded(Obsdisp);

[Obsdisp_mean, Obsdisp_cilog] = errbar4shadedlog(Obsdisp);

[Obsdispzon_mean, Obsdispzon_cilog] = errbar4shadedlog(Obsdisp_zon);
[Obsdispmer_mean, Obsdispmer_cilog] = errbar4shadedlog(Obsdisp_mer);

%% Trial figures

figure(1)

subplot(121)
for i =1:2:length(distance_class)
    for j =1
        loglog(T, Moddisp(:,j,i))
        hold all
        loglog(T, Obsdisp_mean(:,j,i), '--')
    end
end

subplot(122)
for i =1:2:length(distance_class)
    for j =2
        loglog(T, Moddisp(:,j,i))
        hold all
        loglog(T, Obsdisp_mean(:,j,i), '--')
    end
end
loglog(T, 10^5*T.^3)
loglog(T, 10^9*T)

%%
close all 
figure(2)

%%subplot(121)
for i =1:2:length(distance_class)
    for j =2
        loglog(T, Moddisp_zon(:,j,i), 'linewidth',2)
        hold all
        loglog(T, Moddisp_mer(:,j,i), 'linewidth',2)
        %loglog(T, Obsdispzon_mean(:,j,i), '--')
        
    end
end

loglog(T, 10^5*T.^3)
loglog(T, 10^9*T)

%%
figure
%%subplot(121)
for i =1:2:length(distance_class)
    for j =2
        loglog(T, Obsdispzon_mean(:,j,i))
        hold all
        loglog(T, Obsdispmer_mean(:,j,i))
        %loglog(T, Obsdispzon_mean(:,j,i), '--')
    end
end

%% 
colors = get(gca,'ColorOrder');

%%
%%%%%%%%%%%%%%%%%%%%%
%%% log-log plots %%%
%%%%%%%%%%%%%%%%%%%%%
%%% Raw plots *******
%% Figure at first depth 
close all 
clear h g
figure('rend','painters','pos',[10 10 800 600])
i=1
for k =1:2:6
    
    h(i)= shadedErrorBar_log(T+0.001, Obsdisp_mean(:,1,i), Obsdisp_cilog(:,:,1,i), ...
                    {'-','linewidth',3,'color',colors(i*2,:)},1)
                hold all
    g(i) = loglog(T+0.001, Moddisp(:,1,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end

Tax= T(5:50);
t3 = loglog(Tax, 10^7*Tax.^3,'--', 'color', [0.5 0.5 0.5]);
t2 = loglog(Tax, 10^8.15*Tax.^2,'-.', 'color', [0.5 0.5 0.5]);
t1 = loglog(Tax, 10^9.3*Tax, 'color', [0.5 0.5 0.5]);


A = legend([h(1).mainLine,h(2).mainLine, h(3).mainLine, g, t3, t2, t1], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km','~t^3','~t^2','~t'}); 
axis([1 100 10^8 3*10^11])
legend boxoff
set(A, 'location','best')
set(gca,'FontSize', 24)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$D^2(t, D_0) (m^2)$', 'Interpreter','Latex')
saveas(gcf,'../figures/rel_disp_shallow.eps', 'epsc')

%% Figure at second depth 
clear h g
figure('rend','painters','pos',[10 10 800 600])
i=1
for k =1:2:6
    
    h(i)= shadedErrorBar_log(T+0.001, Obsdisp_mean(:,2,i), Obsdisp_cilog(:,:,2,i), ...
                    {'-','linewidth',3,'color',colors(i*2,:)},1)
                hold all
    g(i) = loglog(T+0.001, Moddisp(:,2,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end

Tax= T(5:50);
loglog(Tax, 10^7*Tax.^3,'--', 'color', [0.5 0.5 0.5])
loglog(Tax, 10^8.15*Tax.^2,'-.', 'color', [0.5 0.5 0.5]);
loglog(Tax, 10^9.3*Tax, 'color', [0.5 0.5 0.5])


%legend([h(1).mainLine,h(2).mainLine, h(3).mainLine, g], ...
%    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) ;
axis([1 100 10^8 3*10^11])
%legend boxoff
%set(A, 'location','best')

set(gca,'FontSize', 24)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$D^2(t, D_0) (m^2)$', 'Interpreter','Latex')
saveas(gcf,'../figures/rel_disp_deep.eps', 'epsc')

%% Raw plots with zon and meridional separated 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Zonal - Meridional Plots %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure at first depth Model
close all  
clear h g
figure('rend','painters','pos',[10 10 800 600])
i=1
for k =1:2:6
    
    g(i) = loglog(T+0.001, Moddisp_zon(:,1,i), '-', 'color', colors(i*2,:), 'linewidth', 3)
    hold all
    loglog(T+0.001, Moddisp_mer(:,1,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end

Tax= T(5:50);
t3 = loglog(Tax, 10^7*Tax.^3,'--', 'color', [0.5 0.5 0.5]);
t2 = loglog(Tax, 10^8*Tax.^2,'-.', 'color', [0.5 0.5 0.5]);
t1 = loglog(Tax, 10^9*Tax, 'color', [0.5 0.5 0.5]);

A=legend({'Zon 11km', 'Mer 11km','Zon 33km', 'Mer 33km','Zon 50km', 'Mer 50km','~t^3','~t^2','~t'})

axis([1 100 10^7 10^11])
legend boxoff
set(A, 'location','best')
set(gca,'FontSize', 24)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$D^2(t, D_0) (m^2)$', 'Interpreter','Latex')
saveas(gcf,'../figures/rel_disp_shallow_zm_mod.eps', 'epsc')

%% Figure at second depth Model
%close all 
clear h g
figure('rend','painters','pos',[10 10 800 600])
i=1
for k =1:2:6
    
    g(i) = loglog(T+0.001, Moddisp_zon(:,2,i), '-', 'color', colors(i*2,:), 'linewidth', 3)
    hold all
    loglog(T+0.001, Moddisp_mer(:,2,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end

Tax= T(5:50);
t3 = loglog(Tax, 10^7*Tax.^3,'--', 'color', [0.5 0.5 0.5]);
t2 = loglog(Tax, 10^8*Tax.^2,'-.', 'color', [0.5 0.5 0.5]);
t1 = loglog(Tax, 10^9*Tax, 'color', [0.5 0.5 0.5]);

%A=legend({'Zon 11km', 'Mer 11km','Zon 33km', 'Mer 33km','Zon 50km', 'Mer 50km','~t^3','~t^2','~t'})

axis([1 100 10^7 10^11])
%legend boxoff
%set(A, 'location','best')
set(gca,'FontSize', 24)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$D^2(t, D_0) (m^2)$', 'Interpreter','Latex')
saveas(gcf,'../figures/rel_disp_deep_zm_mod.eps', 'epsc')

%% Figure at first depth Model
close all  
clear h g
figure('rend','painters','pos',[10 10 800 600])

i=1
for k =1:2:6
    
    g(i) = loglog(T+0.001, Obsdispzon_mean(:,1,i), '-', 'color', colors(i*2,:), 'linewidth', 3)
    hold all
    loglog(T+0.001, Obsdispmer_mean(:,1,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end

Tax= T(5:50);
t3 = loglog(Tax, 10^7*Tax.^3,'--', 'color', [0.5 0.5 0.5]);
t2 = loglog(Tax, 10^8*Tax.^2,'-.', 'color', [0.5 0.5 0.5]);
t1 = loglog(Tax, 10^9*Tax, 'color', [0.5 0.5 0.5]);

A=legend({'Zon 10-15km', 'Mer 10-15km','Zon 30-35km', 'Mer 30-35km','Zon 50-55km', 'Mer 50-55km','~t^3','~t^2','~t'})

axis([1 100 10^7.8 10^11])
legend boxoff
set(A, 'location','best')
set(gca,'FontSize', 24)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$D^2(t, D_0) (m^2)$', 'Interpreter','Latex')
saveas(gcf,'../figures/rel_disp_shallow_zm_obs.eps', 'epsc')

%% Figure at second depth Model
close all  
clear h g
figure('rend','painters','pos',[10 10 800 600])

i=1
for k =1:2:6
    
    g(i) = loglog(T+0.001, Obsdispzon_mean(:,2,i), '-', 'color', colors(i*2,:), 'linewidth', 3)
    hold all
    loglog(T+0.001, Obsdispmer_mean(:,2,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end

Tax= T(5:50);
t3 = loglog(Tax, 10^7*Tax.^3,'--', 'color', [0.5 0.5 0.5]);
t2 = loglog(Tax, 10^8*Tax.^2,'-.', 'color', [0.5 0.5 0.5]);
t1 = loglog(Tax, 10^9*Tax, 'color', [0.5 0.5 0.5]);


axis([1 100 10^7.8 10^11])

set(gca,'FontSize', 24)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$D^2(t, D_0) (m^2)$', 'Interpreter','Latex')
saveas(gcf,'../figures/rel_disp_deep_zm_obs.eps', 'epsc')


