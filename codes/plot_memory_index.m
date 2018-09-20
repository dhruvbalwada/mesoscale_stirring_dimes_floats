% Dhruv Balwada
% 3 September 2018
% plot_memory_index.m 
% 
% Memory Index (position Correlation)
% Equation 5.1 in Foussard et al 2017 
% M(t) = <Y.Y_0>/ Y_0<Y^2>^0.5
% Where Y is the separation vector, and Y_o is the initial separation
% vector


%% model memory

ndays=100; 
Mmodel1 = memory_index(model_sep(1).sep_final, ndays); 
Mmodel2 = memory_index(model_sep(2).sep_final, ndays); 

%% observations memory 

Mobs1 = memory_index(obs_sep(1).sep, 100);
Mobs2 = memory_index(obs_sep(2).sep, 100);

for i = 1:500
    y = randsample(pairs_ncon(1),pairs_ncon(1),true);
    
    Mobs1bs(1:ndays,i) = memory_index(obs_sep(1).sep(y), 100);
end


for i = 1:500
    y = randsample(pairs_ncon(2),pairs_ncon(2),true);
    
    Mobs2bs(1:ndays,i) = memory_index(obs_sep(2).sep(y), 100);
end

for i =1:ndays
    Mobs1ci(:,i) = prctile(Mobs1bs(i,:), [95, 5]);
    Mobs2ci(:,i) = prctile(Mobs2bs(i,:), [95, 5]);
end

%%
Mobs1cir(1,:) = Mobs1ci(1,:)' - Mobs1; 
Mobs1cir(2,:) = -Mobs1ci(2,:)' + Mobs1;
Mobs2cir(1,:) = Mobs2ci(1,:)' - Mobs2; 
Mobs2cir(2,:) = -Mobs2ci(2,:)' + Mobs2;

%%
colors =get(gca,'colororder');

%% 

close all
figure
o(1)=shadedErrorBar_semilogx(T+0.001, Mobs1, Mobs1cir, {'-','linewidth',3,'color',colors(1,:)},1);
hold all
o(2)=shadedErrorBar_semilogx(T+0.001, Mobs2, Mobs2cir, {'-','linewidth',3,'color',colors(2,:)},1);
m(1) = semilogx(T+0.001,Mmodel1, '--','linewidth', 3, 'color',colors(3,:));
m(2) = semilogx(T+0.001,Mmodel2, '--', 'linewidth', 3, 'color',colors(4,:));

legend([o(1).mainLine, o(2).mainLine, m], {'500-1000m', '1000-1800m', '750m', '1500m'})
axis([1 100 -0.1 1])

set(gca,'fontsize',20)
xlabel('$t$ (Days)', 'Interpreter','Latex') 

ylabel('$<\bf{y}$ $(t)$.$\bf{y_o}>/$ $y_o<y^2>^{1/2}$', 'Interpreter','Latex')

saveas(gca,'../figures/memory.pdf')

