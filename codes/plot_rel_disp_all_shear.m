%% Code in progress .... 
% Initial tests showed that the system was not like shear dispersion 
% the zonal component did not show a marked T^3 growth. So we gave up
% trying further. 

figure 

for i=1:2:6
for j=2
loglog(T, (Moddisp_zon(:,j,i) - Moddisp_zon(1,j,i)).^(1/3))
hold all
end
end
loglog(T, 100*T,'-','color',[0.5 0.5 0.5], 'linewidth',1)
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