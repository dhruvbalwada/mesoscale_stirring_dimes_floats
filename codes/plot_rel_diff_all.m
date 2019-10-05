% Dhruv Balwada
% 20 September 2018

% calculate derivative of the mean 
for i=1:length(distance_class)
    for j=1:2
        
        Mod_diff(:,j,i) = rel_diff(Moddisp(:,j,i), 100);
        %Obs_diff(:,j,i) = rel_diff(nanmoving_average(Obsdisp_mean(:,j,i),1), 100);
    end
end

for i=1:length(distance_class)
    for j=1:2
        for n =1:500
           Obs_diffbs(:,n,j,i) = rel_diff(nanmoving_average(Obsdisp(:,n,j,i),1), 100);
        end
     end
end

%% calculate mean of the derivative

% model
Mod_diff2 = nan*Mod_diff; 
for i =1:length(distance_class) 
    for j =1:2
        
            Mod_diff2(2:end-1,j,i) = nanmean((Mod_disp(j,i).disp(3:end,:) - Mod_disp(j,i).disp(1:end-2,:))/24/3600/2/2,2);
        
    end
end

%% obs using the raw dispersion
Obs_diff2 = nan*Mod_diff; 

for i=1:length(distance_class)
    for j=1:2
        
        Obs_diff2(2:end-1,j,i) = nanmean((Obs_disp(j,i).disp(3:end,:) - Obs_disp(j,i).disp(1:end-2,:))/24/3600/2/2,2);
        %Obs_diff(:,j,i) = rel_diff(nanmoving_average(Obsdisp_mean(:,j,i),1), 100);
    end
end

%% obs using the bootstrapped dispersion 

Obs_diff2bs = nan*ones(100,500,2,6);
smdays =3;
for i=1:length(distance_class)
    for j=1:2
        for n=1:500 
            Obs_diff2bs(2:end-1,n,j,i) =  ...
                nanmean((nanmoving_average(Obsbs(n,j,i).disp(3:end,:), smdays,1) ...
                - nanmoving_average(Obsbs(n,j,i).disp(1:end-2,:),smdays,1))/24/3600/2/2,2);
        end
    end
end 
%%
[Obsdiff_mean, Obsdiff_ci] = errbar4shaded(Obs_diffbs);
[Obsdiff2_mean, Obsdiff2_ci] = errbar4shadedlog(Obs_diff2bs);

%% Figure rel diff 
% figure shallow 
figure('rend','painters','pos',[10 10 800 600])

i=1;
for k=1:2:length(distance_class)
    for j=1

        h(i) = shadedErrorBar_log(nanmoving_average(Obsdisp_mean(:,j,k).^0.5/1e3,smdays,1), ...
            Obsdiff2_mean(:,j,k), Obsdiff2_ci(:,:,j,k), ...
            {'-','linewidth',3,'color',colors(i*2,:)},1);
        hold all
        g(i) = loglog(Moddisp(:,j,k).^0.5/1e3, Mod_diff(:,j,k),'--','color', colors(i*2,:), 'linewidth', 3);
        i=i+1;
    end
end
axis([10 400 90 30e3])
rax = [10:300] ;

r2 = loglog(rax, 3*rax.^2, 'color',[0.5 0.5 0.5]);
r43 = loglog(rax, 30*rax.^(4/3), '--', 'color', [0.5 0.5 0.5]);


A = legend([h(1).mainLine,h(2).mainLine, h(3).mainLine, g, r2, r43], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km','~r^2','~r^{4/3}'}); 
legend boxoff
set(A, 'location','best')
set(gca,'FontSize', 24)
xlabel('$D$ (km)', 'Interpreter','Latex') 
ylabel('$\kappa(D) (m^2/s)$', 'Interpreter','Latex')
saveas(gcf,'../figures/rel_diff_shallow.eps', 'epsc')
%% Figure rel diff 
% figure deep 
figure('rend','painters','pos',[10 10 800 600])

i=1;
for k=1:2:length(distance_class)
    for j=2

        h(i) = shadedErrorBar_log(nanmoving_average(Obsdisp_mean(:,j,k).^0.5/1e3,smdays,1), ...
            Obsdiff2_mean(:,j,k), Obsdiff2_ci(:,:,j,k), ...
            {'-','linewidth',3,'color',colors(i*2,:)},1);
        hold all
        g(i) = loglog(Moddisp(:,j,k).^0.5/1e3, Mod_diff(:,j,k),'--','color', colors(i*2,:), 'linewidth', 3);
        i=i+1;
    end
end
axis([10 400 20 20e3])
rax = [10:300] ;

r2 = loglog(rax, 3*rax.^2, 'color',[0.5 0.5 0.5]);
r43 = loglog(rax, 30*rax.^(4/3), '--', 'color', [0.5 0.5 0.5]);


%A = legend([h(1).mainLine,h(2).mainLine, h(3).mainLine, g, r2, r43], ...
%    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km','~r^2','~r^{4/3}'}); 
%legend boxoff
%set(A, 'location','best')
set(gca,'FontSize', 24)
xlabel('$D$ (km)', 'Interpreter','Latex') 
ylabel('$\kappa(D) (m^2/s)$', 'Interpreter','Latex')
saveas(gcf,'../figures/rel_diff_deep.eps', 'epsc')

%% Char time scale - plots vs D 
% first depth
clear h g
figure('rend','painters','pos',[10 10 800 600])
i=1;
for k=1:2:length(distance_class)
    for j=1
        
        h(i) = semilogx(nanmoving_average(Obsdisp_mean(:,j,k).^0.5/1e3,smdays,1), ...
            nanmoving_average(Obsdisp_mean(:,j,k),smdays,1)./Obsdiff2_mean(:,j,k)/24/3600,...
            'o-','color', colors(i*2,:), 'linewidth', 3);
         hold all
        g(i) = semilogx(Moddisp(:,j,k).^0.5/1e3, Moddisp(:,j,k)./Mod_diff(:,j,k)/24/3600,...
            '--','color', colors(i*2,:), 'linewidth', 3);
        %semilogx(T, Moddisp(:,j,i)./Mod_diff(:,j,i)/24/3600,'.-')
       
        i = i+1;
    end
end
axis([10 400 0 150])
rax = 0:400; 
semilogx(rax, 8*rax./rax, '-', 'linewidth',1, 'color',[0.5 0.5 0.5])
A = legend([h, g], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}); 
legend boxoff
set(A, 'location','best')
set(gca,'FontSize', 24)
xlabel('$D$ (km)', 'Interpreter','Latex') 
ylabel('$\tau (Days)$', 'Interpreter','Latex')
saveas(gcf,'../figures/char_time_shallow.eps', 'epsc')

%% second depth 
clear h g
figure('rend','painters','pos',[10 10 800 600])
i=1;
for k=1:2:length(distance_class)
    for j=2
        
        h(i) = semilogx(nanmoving_average(Obsdisp_mean(:,j,k).^0.5/1e3,smdays,1), ...
            nanmoving_average(Obsdisp_mean(:,j,k),smdays,1)./Obsdiff2_mean(:,j,k)/24/3600,...
            'o-','color', colors(i*2,:), 'linewidth', 3);
         hold all
        g(i) = semilogx(Moddisp(:,j,k).^0.5/1e3, Moddisp(:,j,k)./Mod_diff(:,j,k)/24/3600,...
            '--','color', colors(i*2,:), 'linewidth', 3);
        %semilogx(T, Moddisp(:,j,i)./Mod_diff(:,j,i)/24/3600,'.-')
       
        i = i+1;
    end
end
axis([10 400 0 150])
semilogx(rax, 15*rax./rax, '-', 'linewidth',1, 'color',[0.5 0.5 0.5])
%A = legend([h, g], ...
%    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'}); 
%legend boxoff
%set(A, 'location','best')
set(gca,'FontSize', 24)
xlabel('$D$ (km)', 'Interpreter','Latex') 
ylabel('$\tau (Days)$', 'Interpreter','Latex')
saveas(gcf,'../figures/char_time_deep.eps', 'epsc')
%% Char time scale  - plots vs T 
figure
for k=1:2:length(distance_class)
    for j=1
        
        semilogx(T, Moddisp(:,j,k)./Mod_diff(:,j,k)/24/3600,'--')
        %semilogx(T, Moddisp(:,j,i)./Mod_diff(:,j,i)/24/3600,'.-')
        hold all
        semilogx(T, ...
            nanmoving_average(Obsdisp_mean(:,j,k),smdays,1)./Obsdiff2_mean(:,j,k)/24/3600)
        
    end
end
axis([1 100 0 150])
