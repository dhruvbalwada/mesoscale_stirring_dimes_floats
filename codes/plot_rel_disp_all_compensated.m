%% Semilog Y axis plots 
% Detect Lin's law 
% to detect if D^2 ~ D_o^2 exp(lambda . t)

%% Figure at first depth 
close all 
clear h g
figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    
  %  h(i)= shadedErrorBar_semilogy(T+0.001, Obsdisp_mean(:,1,i), Obsdisp_cilog(:,:,1,i), ...
  %                  {'-','linewidth',3,'color',colors(i*2,:)},1)
     h(i)= semilogy(T, Obsdisp_mean(:,1,i), ...
                   'o-','linewidth',3,'color',colors(i*2,:))
   
                hold all
    g(i) = semilogy(T, Moddisp(:,1,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end
A = legend([h(1),h(2), h(3), g], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) 
axis([1 100 10^8 3*10^11])
legend boxoff
set(A, 'location','best')
set(gca,'FontSize', 24)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$D^2(t, D_0)$', 'Interpreter','Latex')
saveas(gcf,'../figures/rel_disp_shallow_exp.eps', 'epsc')

%% Figure at second depth 
%close all 
clear h g
figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    
  %  h(i)= shadedErrorBar_semilogy(T+0.001, Obsdisp_mean(:,1,i), Obsdisp_cilog(:,:,1,i), ...
  %                  {'-','linewidth',3,'color',colors(i*2,:)},1)
     h(i)= semilogy(T, Obsdisp_mean(:,2,i), ...
                   'o-','linewidth',3,'color',colors(i*2,:))
   
                hold all
    g(i) = semilogy(T, Moddisp(:,2,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end
A = legend([h(1),h(2), h(3), g], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) 
axis([1 100 10^8 3*10^11])
legend boxoff
set(A, 'location','best')
set(gca,'FontSize', 24)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$D^2(t, D_0)$', 'Interpreter','Latex')

saveas(gcf,'../figures/rel_disp_deep_exp.eps', 'epsc')

%% Try and detect ballistic growth 
% D^2 = Do^2 + ct^2 
% log(D^2 - Do^2) = 2log(ct)

%% Figure at first depth (no errorbars)
%close all 
clear h g
figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    
    h(i)= loglog(T, Obsdisp_mean(:,1,i) - Obsdisp_mean(1,1,i), ...
                    'o-','linewidth',3,'color',colors(i*2,:))
                hold all
    g(i) = loglog(T, Moddisp(:,1,i) - Moddisp(1,1,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end

t2 = loglog(T, 10^8*T.^2,'--','color',[0.5 0.5 0.5], 'linewidth',1)

A = legend([h(1),h(2), h(3), g, t2], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km', '~t^2'}) ;
axis([1 100 10^6.9 10^11])
legend boxoff
set(A, 'location','northwest')
set(gca,'FontSize', 24)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$[D^2(t, D_0) - D_0^2] (m^2)$', 'Interpreter','Latex')
saveas(gcf,'../figures/rel_disp_shallow_bal.eps', 'epsc')

%% Figure at second depth (no errorbars)
close all 
clear h g
figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    
    h(i)= loglog(T, Obsdisp_mean(:,2,i) - Obsdisp_mean(1,2,i), ...
                    'o-','linewidth',3,'color',colors(i*2,:))
                hold all
    g(i) = loglog(T, Moddisp(:,2,i) - Moddisp(1,2,i), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end

t2 = loglog(T, 10^7.5*T.^2,'--','color',[0.5 0.5 0.5], 'linewidth',1)

%A = legend([h(1),h(2), h(3), g, t2], ...
%    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km', '~t^2'}) ;
axis([1 100 10^6.9 10^11])
%legend boxoff
%set(A, 'location','best')
set(gca,'FontSize', 24)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$[D^2(t, D_0) - D_0^2] (m^2)$', 'Interpreter','Latex')
saveas(gcf,'../figures/rel_disp_deep_bal.eps', 'epsc')


%% compensated plots first depth
close all 
clear h g
figure('rend','painters','pos',[10 10 400 300])
i=1;

for k =1:2:6
    
    h(i)= semilogx(T, (Obsdisp_mean(:,1,i) - Obsdisp_mean(1,1,i))./T'.^2, ...
                    'o-','linewidth',3,'color',colors(i*2,:));
                hold all
    g(i) = semilogx(T, (Moddisp(:,1,i) - Moddisp(1,1,i))./T'.^2, '--', 'color', colors(i*2,:), 'linewidth', 3);
    
    i=i+1;
end

%A = legend([h, g], ...
%    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) ;
axis([1 100 0 9*10^7])
%legend boxoff
%set(A, 'location','best')
set(gca,'FontSize', 18)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$[D^2(t, D_0) - D_0^2]/t^2$', 'Interpreter','Latex')
saveas(gcf,'../figures/compen_rel_disp_shallow_bal.eps', 'epsc')

%% compensated plots second depth
close all 
clear h g
figure('rend','painters','pos',[10 10 400 300])
i=1;

for k =1:2:6
    
    h(i)= semilogx(T, (Obsdisp_mean(:,2,i) - Obsdisp_mean(1,2,i))./T'.^2, ...
                    '-','linewidth',3,'color',colors(i*2,:));
                hold all
    g(i) = semilogx(T, (Moddisp(:,2,i) - Moddisp(1,2,i))./T'.^2, '--', 'color', colors(i*2,:), 'linewidth', 3);
    
    i=i+1;
end

%A = legend([h, g], ...
%    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) ;
axis([1 100 0 3*10^7])
%legend boxoff
%set(A, 'location','best')
set(gca,'FontSize', 18)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$[D^2(t, D_0) - D_0^2]/t^2$', 'Interpreter','Latex')
saveas(gcf,'../figures/compen_rel_disp_deep_bal.eps', 'epsc')

%% Check Richardson's law (cubic growth). 

%% Figure at first depth (no errorbars)
close all 
clear h g
figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    
    h(i)= loglog(T, Obsdisp_mean(:,1,i).^(1/3) - Obsdisp_mean(1,1,i).^(1/3), ...
                    'o-','linewidth',3,'color',colors(i*2,:))
                hold all
    g(i) = loglog(T, Moddisp(:,1,i).^(1/3) - Moddisp(1,1,i).^(1/3), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end

t1= loglog(T, 100*T,'-','color',[0.5 0.5 0.5], 'linewidth',1)
axis([1 100 3 5*10^3])
A = legend([h(1),h(2), h(3), g, t1], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km', '~t^1'}) 
%axis([1 100 10^6.9 10^11])
legend boxoff
set(A, 'location','northwest')
set(gca,'FontSize', 24)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('[$D^{2/3}(t, D_0) - D_0^{2/3}] (m^{2/3})$', 'Interpreter','Latex')

saveas(gcf,'../figures/rel_disp_shal_rich.eps', 'epsc')

%% Figure at second depth (no errorbars)
%close all 
clear h g
figure('rend','painters','pos',[10 10 900 600])
i=1
for k =1:2:6
    
    h(i)= loglog(T, Obsdisp_mean(:,2,i).^(1/3) - Obsdisp_mean(1,2,i).^(1/3), ...
                    'o-','linewidth',3,'color',colors(i*2,:))
                hold all
    g(i) = loglog(T, Moddisp(:,2,i).^(1/3) - Moddisp(1,2,i).^(1/3), '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end

loglog(T, 100*T,'-','color',[0.5 0.5 0.5], 'linewidth',1)

%A = legend([h(1),h(2), h(3), g], ...
%    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}) 
axis([1 100 3 5*10^3])
%legend boxoff
%set(A, 'location','best')
set(gca,'FontSize', 24)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('[$D^{2/3}(t, D_0) - D_0^{2/3}] (m^{2/3})$', 'Interpreter','Latex')
saveas(gcf,'../figures/rel_disp_deep_rich.eps', 'epsc')

%% compensated plot depth 1 

close all 
clear h g
figure('rend','painters','pos',[10 10 400 300])
i=1
depth=1
for k =1:2:6
    
    h(i)= semilogx(T, (Obsdisp_mean(:,depth,i).^(1/3) - Obsdisp_mean(1,depth,i).^(1/3))./T', ...
                    'o-','linewidth',3,'color',colors(i*2,:))
                hold all
    g(i) = semilogx(T, (Moddisp(:,depth,i).^(1/3) - Moddisp(1,depth,i).^(1/3))./T', '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end

axis([1 100 0 110])

set(gca,'FontSize', 18)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$(D^{2/3}(t, D_0) - D_0^{2/3})/t $', 'Interpreter','Latex')
saveas(gcf,'../figures/compen_rel_disp_shallow_rich.eps', 'epsc')

%% compensated plot depth 2

close all 
clear h g
figure('rend','painters','pos',[10 10 400 300])
i=1
depth=2
for k =1:2:6
    
    h(i)= semilogx(T, (Obsdisp_mean(:,depth,i).^(1/3) - Obsdisp_mean(1,depth,i).^(1/3))./T', ...
                    'o-','linewidth',3,'color',colors(i*2,:))
                hold all
    g(i) = semilogx(T, (Moddisp(:,depth,i).^(1/3) - Moddisp(1,depth,i).^(1/3))./T', '--', 'color', colors(i*2,:), 'linewidth', 3)
    
    i=i+1
end

axis([1 100 0 60])


set(gca,'FontSize', 18)
xlabel('$t$ (Days)', 'Interpreter','Latex') 
ylabel('$(D^{2/3}(t, D_0) - D_0^{2/3})/t $', 'Interpreter','Latex')

saveas(gcf,'../figures/compen_rel_disp_deep_rich.eps', 'epsc')