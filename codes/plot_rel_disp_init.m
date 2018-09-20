% Dhruv Balwada
% 3 September 2018
% plot_rel_disp_init.m 
% 
% In the small time regime we expect <y^2> = <y_o^2> + <dv_o.dv_o>t^2
% this is referred to as the ballistic regime 

%% model

Mod_disp1 = rel_disp(model_sep(1).sep_final,ndays,1);
Mod_disp2 = rel_disp(model_sep(2).sep_final,ndays,1);

%% 
Mod_dispnorm1 = Mod_disp1.avdisp - Mod_disp1.avdisp(1); 
Mod_dispnorm2 = Mod_disp2.avdisp - Mod_disp2.avdisp(1); 

%% Observations

Obs_disp1 = rel_disp(obs_sep(1).sep, ndays,1);
Obs_disp2 = rel_disp(obs_sep(2).sep, ndays,1);

Obs_dispnorm1 = Obs_disp1.avdisp - Obs_disp1.avdisp(1); 
Obs_dispnorm2 = Obs_disp2.avdisp - Obs_disp2.avdisp(1); 

%%
for i = 1:500
    y = randsample(pairs_ncon(1),pairs_ncon(1),true);
    
    Obs1bs = rel_disp(obs_sep(1).sep(y), ndays,1);
    
    Obs_dispnorm1bs(1:ndays, i) = Obs1bs.avdisp - Obs1bs.avdisp(1); 
end

for i = 1:500
    y = randsample(pairs_ncon(2),pairs_ncon(2),true);
    
    Obs2bs = rel_disp(obs_sep(2).sep(y), ndays,1);
    
    Obs_dispnorm2bs(1:ndays, i) = Obs2bs.avdisp - Obs2bs.avdisp(1); 
end

%%

for i =1:ndays
    Obs_dispnorm1ci(:,i) = prctile(Obs_dispnorm1bs(i,:), [95, 5]);
    Obs_dispnorm2ci(:,i) = prctile(Obs_dispnorm2bs(i,:), [95, 5]);
end

%%
Obs_disp1ci(1,:) = 0.434*(Obs_dispnorm1ci(1,:)' - Obs_dispnorm1) ...
                ./Obs_dispnorm1; 
Obs_disp1ci(2,:) = 0.434*(-Obs_dispnorm1ci(2,:)' + Obs_dispnorm1) ...
                ./Obs_dispnorm1; 
            
Obs_disp2ci(1,:) = 0.434*(Obs_dispnorm2ci(1,:)' - Obs_dispnorm2) ...
                ./Obs_dispnorm2; 
Obs_disp2ci(2,:) = 0.434*(-Obs_dispnorm2ci(2,:)' + Obs_dispnorm2) ...
                ./Obs_dispnorm2; 


%%
colors =get(gca,'colororder');
            
%% error bar
close all
figure
o(1)=shadedErrorBar_log(T(2:end), Obs_dispnorm1(2:end), Obs_disp1ci(:,2:end), {'o-','linewidth',3,'color',colors(1,:)},1);
hold all
o(2)=shadedErrorBar_log(T(2:end), Obs_dispnorm2(2:end), Obs_disp2ci(:,2:end), {'-','linewidth',3,'color',colors(2,:)},1);

m(1) = loglog(T, Mod_dispnorm1, '--', 'linewidth',3,'color',colors(3,:));
m(2) = loglog(T, Mod_dispnorm2, '--', 'linewidth',3,'color',colors(4,:));

th(1) = loglog(T, 10^(7.6)*T.^2, '--', 'linewidth',2,'color',[0.5 0.5 0.5]);

set(gca, 'fontsize', 20) 
axis([1 100 10^6.7 10^11.5])

legend([o(1).mainLine, o(2).mainLine, m, th], {'500-1000m', '1000-1800m', '750m', '1500m','~t^{2}'},'location','best')

xlabel('$t$ (Days)', 'Interpreter', 'Latex')
ylabel('$<y^2>(t) - y_o^2$', 'Interpreter', 'Latex')

saveas(gca,'../figures/ballistic_disp.pdf')

%%
close all
figure
loglog(T, Mod_dispnorm1,'linewidth',2) 
hold all 
loglog(T, Mod_dispnorm2,'linewidth',2) 
loglog(T, Obs_dispnorm1,'.-','linewidth',2) 
loglog(T, Obs_dispnorm2,'.-','linewidth',2) 
%loglog(T+0.001, 10^(7.3)*T.^2, '--', 'linewidth',2,'color',[0.5 0.5 0.5])
