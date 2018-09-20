% Dhruv Balwada
% 3 September 2018
% plot_chance_correlation.m 
% 
% Memory Index (position Correlation)
% Equation 5.1 in Foussard et al 2017 
% M(t) = <Y.Y_0>/ Y_0<Y^2>^0.5
% Where Y is the separation vector, and Y_o is the initial separation
% vector


%% model memory

ndays=100; 
Cmodel1 = chance_correlation(model_sep(1).sep_final, ndays); 
Cmodel2 = chance_correlation(model_sep(2).sep_final, ndays); 

%% observations memory 

Cobs1 = chance_correlation(obs_sep(1).sep, 100);
Cobs2 = chance_correlation(obs_sep(2).sep, 100);

for i = 1:500
    y = randsample(pairs_ncon(1),pairs_ncon(1),true);
    
    Cobs1bs(1:ndays,i) = chance_correlation(obs_sep(1).sep(y), 100);
end


for i = 1:500
    y = randsample(pairs_ncon(2),pairs_ncon(2),true);
    
    Cobs2bs(1:ndays,i) = chance_correlation(obs_sep(2).sep(y), 100);
end

for i =1:ndays
    Cobs1ci(:,i) = prctile(Cobs1bs(i,:), [95, 5]);
    Cobs2ci(:,i) = prctile(Cobs2bs(i,:), [95, 5]);
end

%%
Cobs1cir(1,:) = Cobs1ci(1,:)' - Cobs1; 
Cobs1cir(2,:) = -Cobs1ci(2,:)' + Cobs1;
Cobs2cir(1,:) = Cobs2ci(1,:)' - Cobs2; 
Cobs2cir(2,:) = -Cobs2ci(2,:)' + Cobs2;

%%
colors =get(gca,'colororder');

%% 

close all
figure
o(1)=shadedErrorBar_semilogx(T+0.001, Cobs1, Cobs1cir, {'-','linewidth',3,'color',colors(1,:)},1);
hold all
o(2)=shadedErrorBar_semilogx(T+0.001, Cobs2, Cobs2cir, {'-','linewidth',3,'color',colors(2,:)},1);
m(1) = semilogx(T+0.001,Cmodel1, '--','linewidth', 3, 'color',colors(3,:));
m(2) = semilogx(T+0.001,Cmodel2, '--', 'linewidth', 3, 'color',colors(4,:));

legend([o(1).mainLine, o(2).mainLine, m], {'500-1000m', '1000-1800m', '750m', '1500m'},'location','best')
axis([1 100 -0.5 0.5])

set(gca,'fontsize',20)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$<\delta$ $\bf{v}$ $(t).$  $\bf{y_o}>$', 'Interpreter','Latex')
saveas(gca,'../figures/chance_correlation.pdf')

