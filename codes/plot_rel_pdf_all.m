% Dhruv Balwada
% 21 September 2018

%1. look at the PDFs of relative separation in log form
%2. Plot kurtosis plots 

%% 
y_0 = nanmean(Mod_disp(2,1).disp(1,:).^0.5); 
xlog = logspace(-2,2.5,20); 
xaxislog = 0.5*(xlog(1:end-1) + xlog(2:end));

%% 

std1 = nanstd(log(Mod_disp(2,1).disp'.^0.5/y_0)); 

%%
clear pdf_lung
for i =1:100
h1=histogram(((Mod_disp(2,1).disp(i,:).^0.5)/y_0).^(1/std1(i)), xlog);
pdf_lung(:,i) = h1.Values./max(h1.Values); 

h2=histogram(((Mod_disp(2,1).disp(i,:).^0.5)/y_0), xlog);
pdf_rich(:,i) = h2.Values./max(h2.Values); 

end

%%
close all 
figure
for i=4:5:40
loglog(xaxislog, pdf_lung(:,i),'.-','linewidth',2)
hold all 
end

figure
for i=4:5:40
loglog(xaxislog, pdf_rich(:,i),'.-','linewidth',2)
hold all 
end


%% look at kurtosis

clear kurt_mod kurt_obs

for i =1:length(distance_class)
    for j =1:2        
            kurt_mod(:,j,i) = nanmean(Mod_disp(j,i).disp'.^2, 1)./nanmean(Mod_disp(j,i).disp', 1).^2  ;
            kurt_obs(:,j,i) = nanmean(Obs_disp(j,i).disp'.^2, 1)./nanmean(Obs_disp(j,i).disp', 1).^2  ;
    end
end

%%
for i=1:length(distance_class)
     for j =1:2
         for n=1:500
            kurt_obsbs(:, n, j, i) = nanmean(Obsbs(n,j,i).disp'.^2, 1)./nanmean(Obsbs(n,j,i).disp', 1).^2  ;
            %kurtosis(Obsbs(n,j,i).disp'.^0.5,1); 
         end
     end
end

%%
[kurt_obs_mean, kurt_obs_ci] = errbar4shaded(kurt_obsbs);

[kurt_obs_mean, kurt_obs_cilog] = errbar4shadedlog(kurt_obsbs);



%% 
close all 
figure('rend','painters','pos',[10 10 800 600])
i=1; 
for k=1:2:length(distance_class) 
    for j=1
        
        shadedErrorBar_semilogy(T, kurt_obs_mean(:,j,k), kurt_obs_cilog(:,:,j,k), ...
            {'-','linewidth',3,'color',colors(i*2,:)},1)
        hold all
        semilogy(T, kurt_mod(:,j,k),'--', 'color', colors(i*2,:), 'linewidth',2)
        i=i+1;
        
    end
end

semilogy(T, 5.6*T./T, '--', 'color',[0.5 0.5 0.5], 'linewidth', 1)

semilogy(T, exp(8*T/fit_parms_obs(j).Tlung), '-', 'color',[0.5 0.5 0.5], 'linewidth', 1) 

%A = legend([h(1).mainLine, h(2).mainLine, h(3).mainLine, g], ...
%    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km'});
    

axis([0 100 1 20])

set(gca, 'FontSize', 24) 
xlabel('$t$ (Days)', 'Interpreter', 'Latex')
ylabel('Kurtosis')
saveas(gcf,'../figures/kurt_shallow.eps', 'epsc')

%%
figure('rend','painters','pos',[10 10 800 600])
i=1; 
clear h g
for k=1:2:length(distance_class) 
    for j=2
        
        h(i) = shadedErrorBar_semilogy(T, kurt_obs_mean(:,j,k), kurt_obs_cilog(:,:,j,k), ...
            {'-','linewidth',3,'color',colors(i*2,:)},1);
        hold all
        g(i) = semilogy(T, kurt_mod(:,j,k),'--', 'color', colors(i*2,:), 'linewidth',2);
        i=i+1;
        
    end
end

rich = semilogy(T, 5.6*T./T, '--', 'color',[0.5 0.5 0.5], 'linewidth', 1);
lung = semilogy(T, exp(8*T/fit_parms_obs(j).Tlung), '-', 'color',[0.5 0.5 0.5], 'linewidth', 1) ;

%semilogy(T, exp(T/2)) 

A = legend([h(1).mainLine, h(2).mainLine, h(3).mainLine, g, rich, lung], ...
    {'Obs 10-15km', 'Obs 30-35km', 'Obs 50-55km','Mod 11km', 'Mod 33km', 'Mod 50km', 'Richardson', 'Lungdren'});
    
set(A, 'location', 'best' , 'FontSize', 18) 
axis([0 100 1 20])
legend boxoff

set(gca, 'FontSize', 24) 
xlabel('$t$ (Days)', 'Interpreter', 'Latex')
ylabel('Kurtosis')
saveas(gcf,'../figures/kurt_deep.eps', 'epsc')


%% 
for k=1:length(distance_class)
for j=1:2
    Tscale_mod(:,j,k) = kurt_mod(2:end-1,j,k)./(kurt_mod(3:end,j,k) - kurt_mod(1:end-2,j,k))/2;
    Tscale_obs(:,j,k) = kurt_obs_mean(2:end-1,j,k)./(kurt_obs_mean(3:end,j,k) - kurt_obs_mean(1:end-2,j,k))/2;
end
end

%%
figure
for k=1:2:length(distance_class)
    for j =2
    plot( Tscale_mod(:,j,k),'o-')
    hold all
    plot( Tscale_obs(:,j,k),'--')
    end
end
axis([0 20 0 20])