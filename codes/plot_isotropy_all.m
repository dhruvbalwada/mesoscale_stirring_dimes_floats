%% Isotropy 

for i =1:length(distance_class)
    for j =1:2
        for n=1:500
            Iso_obs_bs(:,n,j,i) = Obsbs(n,j,i).avdispzon./Obsbs(n,j,i).avdispmer; 
            Robs_bs(:,n,j,i) = Obsbs(n,j,i).avdisp.^0.5; 
        end
    end
end

%%
[Iso_obs, Iso_obsci] = errbar4shaded(Iso_obs_bs); 
Robs = squeeze(nanmean(Robs_bs,2)); 

for i =1:length(distance_class)
    for j =1:2
            Iso_mod(:,j,i) = Mod_disp(j,i).avdispzon./Mod_disp(j,i).avdispmer; 
            Rmod(:,j,i) = Mod_disp(j,i).avdisp.^0.5; 
    end
end

%%
colors = get(gca,'ColorOrder');

%% Iso plots with time 
%% Figure at first depth 
close all 
clear h g
figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    
    h(i)= shadedErrorBar_semilogx(T+0.001, Iso_obs(:,1,i), Iso_obsci(:,:,1,i), ...
                    {'-','linewidth',3,'color',colors(i*2,:)},1)
                hold all
    g(i) = semilogx(T+0.001, Iso_mod(:,1,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end
legend([h(1).mainLine,h(2).mainLine, h(3).mainLine, g], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) 
axis([1 100 0 8])

set(gca,'FontSize', 20)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$<D_x^2>/<D_y^2>$', 'Interpreter','Latex')

%% Figure at second depth 
close all 
clear h g
figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    
    h(i)= shadedErrorBar_semilogx(T+0.001, Iso_obs(:,2,i), Iso_obsci(:,:,2,i), ...
                    {'-','linewidth',3,'color',colors(i*2,:)},1)
                hold all
    g(i) = semilogx(T+0.001, Iso_mod(:,2,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end
legend([h(1).mainLine,h(2).mainLine, h(3).mainLine, g], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) 
axis([1 100 0 3])

set(gca,'FontSize', 20)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$<D_x^2>/<D_y^2>$', 'Interpreter','Latex')

%% Isoplots with separation 
%% Figure at first depth 
close all 
clear h g
figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    
    h(i)= shadedErrorBar_semilogx(Robs(:,1,i)/1e3, Iso_obs(:,1,i), Iso_obsci(:,:,1,i), ...
                    {'-','linewidth',3,'color',colors(i*2,:)},1)
                hold all
    g(i) = semilogx(Rmod(:,1,i)/1e3, Iso_mod(:,1,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end
legend([h(1).mainLine,h(2).mainLine, h(3).mainLine, g], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) 
axis([10 300 0 5])

set(gca,'FontSize', 20)
xlabel('$r$ (km)', 'Interpreter','Latex') 
ylabel('$<D_x^2>/<D_y^2>$', 'Interpreter','Latex')

%% Figure at second depth 
close all 
clear h g
figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    
    h(i)= shadedErrorBar_semilogx(Robs(:,2,i)/1e3, Iso_obs(:,2,i), Iso_obsci(:,:,2,i), ...
                    {'-','linewidth',3,'color',colors(i*2,:)},1)
                hold all
    g(i) = semilogx(Rmod(:,2,i)/1e3, Iso_mod(:,2,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end
legend([h(1).mainLine,h(2).mainLine, h(3).mainLine, g], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) 
axis([10 300 0 2])

set(gca,'FontSize', 20)
xlabel('$r$ (km)', 'Interpreter','Latex') 
ylabel('$<D_x^2>/<D_y^2>$', 'Interpreter','Latex')
%%
