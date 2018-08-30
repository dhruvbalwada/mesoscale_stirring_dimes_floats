
function [dispersion] = rel_disp(sep,ndays, figno)

npairs = length(sep);
X = nan*ones(ndays,npairs);
Y = nan*ones(ndays,npairs);
dist = nan*ones(ndays,npairs);

for i = 1:npairs
    len = length(sep(i).X);
    if ndays<=len
        X(1:ndays,i) = sep(i).X(1:ndays);
        Y(1:ndays,i) = sep(i).Y(1:ndays);
        dist(1:ndays,i) = sep(i).dist(1:ndays);
    else
        
    end
end

mer_disp = nanmean(Y.^2,2);
mer_disp_std = nanstd(Y.^2,1,2);
zon_disp = nanmean(X.^2,2);
mer_disp_std = nanstd(X.^2,1,2);
disp = nanmean(dist.^2,2);
disp_std = nanstd(dist.^2,1,2);

dispersion.disp = dist.^2;
dispersion.dipzon = X.^2; 
dispersion.dipmer = Y.^2; 
dispersion.avdispzon = zon_disp;
dispersion.avdispmer = mer_disp;
dispersion.avdisp = disp;

flag = 0;
if flag == 1
    figure(figno-1)
    plot(sort(dist(1,:)),'.','markersize',10)
    hold all
    set(gca,'fontsize',20)
    xlabel('Pair no')
    ylabel('Distance')
    %% log plots
    
    figure(figno)
    loglog([0:ndays-1] , disp./(10^6),'.-','linewidth',2)
    hold all
    set(gca,'fontsize',20)
    xlabel('Time(days)')
    ylabel('Dispersion(m^2)')
    
    figure(figno+1)
    loglog([0:ndays-1] , zon_disp./(10^6),'.-','linewidth',2)
    hold all
    set(gca,'fontsize',20)
    xlabel('Time(days)')
    ylabel('Zonal Dispersion(m^2)')
    
    figure(figno+2)
    loglog([0:ndays-1] , mer_disp./(10^6),'.-','linewidth',2)
    hold all
    set(gca,'fontsize',20)
    xlabel('Time(days)')
    ylabel('Meridional Dispersion(m^2)')
    
    %% semilog y plots
    figure(figno+3)
    semilogy([0:ndays-1] , disp./(10^6),'.-','linewidth',2)
    hold all
    set(gca,'fontsize',20)
    xlabel('Time(days)')
    ylabel('Dispersion(m^2)')
    
    figure(figno+4)
    semilogy([0:ndays-1] , zon_disp./(10^6),'.-','linewidth',2)
    hold all
    set(gca,'fontsize',20)
    xlabel('Time(days)')
    ylabel('Zonal Dispersion(m^2)')
    
    figure(figno+5)
    semilogy([0:ndays-1] , mer_disp./(10^6),'.-','linewidth',2)
    hold all
    set(gca,'fontsize',20)
    xlabel('Time(days)')
    ylabel('Meridional Dispersion(m^2)')
    
    %%
    figure(figno+6)
    plot(zon_disp./mer_disp)
    hold all
    set(gca,'fontsize',20)
    xlabel('Time(days)')
    ylabel('Isotropy(t)')
    
    figure(figno+7)
    plot(disp.^0.5, zon_disp./mer_disp)
    hold all
    set(gca,'fontsize',20)
    xlabel('Separation(m)')
    ylabel('Isotropy(r)')
    
end
%%

