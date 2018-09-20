% Last Modified 29 August 2018
% Dhruv Balwada
% 
% time_dep_obs.m
% 
% This script calculates the time dependent measures from the dimes floats 
% 1) Relative dispersion
% 2) 

clear all
close all

traj = load ('../data/traj_09_03_13_clean.mat');

% delete sections in the Scotia Sea.
traj.Yc(traj.Xc>-82) = NaN;
traj.Xc(traj.Xc>-82) = NaN;

%% calculate the separation time series for all possible pairs.
sep = calculate_seperation_timeseries(traj);

%% Specify pressure criterion, minimum days, and distance classes.
plevel = [500 1000 1800];

ndays  = 100;
pdiff  = 200;

%% Find pairs corresponding to the criterion. 

distance_class_obs = [15,20]*1000;
 
[pairs_ncon(1),pos_ncon(1)] = find_pairs(sep, distance_class_obs, plevel(1:2), pdiff, ndays);
[pairs_ncon(2),pos_ncon(2)] = find_pairs(sep, distance_class_obs, plevel(2:3), pdiff, ndays);

T = 0:ndays-1;

obs_sep(1) = pos_ncon(1);
obs_sep(2) = pos_ncon(2);

%%
dist_class1 = [10, 15]*1000; 
dist_class2 = [20, 25]*1000; 
dist_class3 = [30, 35]*1000;
dist_class4 = [40, 45]*1000; 
dist_class5 = [50, 55]*1000; 
%%
[pairs(1), pos(1)] = find_pairs(sep, dist_class1, plevel(2:3), pdiff, ndays);
[pairs(2), pos(2)] = find_pairs(sep, dist_class2, plevel(2:3), pdiff, ndays);
[pairs(3), pos(3)] = find_pairs(sep, dist_class3, plevel(2:3), pdiff, ndays);
[pairs(4), pos(4)] = find_pairs(sep, dist_class4, plevel(2:3), pdiff, ndays);
[pairs(5), pos(5)] = find_pairs(sep, dist_class5, plevel(2:3), pdiff, ndays);

%% Number of pairs as a function of time 

count = zeros(ndays,2); 

for k =1:2

for i = 1:ndays
    for j=1:pairs_ncon(k)
        if ~isnan(pos_ncon(k).sep(j).dist(i))
            count(i,k) = count(i,k)+1;
        end
    end
end
end

%% 
figure
plot(count, '-', 'linewidth',2)
legend('500-1000m', '1000-1800m', 'location', 'best')
xlabel('Days')
ylabel('No. Pairs')
set(gca,'fontsize',20)

saveas(gca,'../figures/no_pairs.pdf')


%% Initial distribution of pairs

for k =1:2 
for i = 1:pairs_ncon 
   ini_dist(k).dist(i) = pos_ncon(k).sep(i).dist(1); 
end
end
    
%%
xdata=(15:20);
figure
[h1] =histogram(ini_dist(1).dist/1e3, xdata);
hold all
[h2] =histogram(ini_dist(2).dist/1e3, xdata);
legend('500-1000m', '1000-1800m', 'location', 'best')
xlabel('Initial Separation (km)')
ylabel('No. Pairs')
set(gca,'fontsize',20)
saveas(gca,'../figures/ini_hist.pdf')


%% 
for i =1:5
    [M(:,i), Mvel(:,i)] = memory_index(pos(i).sep, 100); 
end

%
figure
for i=1:5
semilogx(Mvel(:,i)/Mvel(1,i),'.-')
hold all
%semilogx(M(:,i),'--')
end

%%
for i =1:5
    C(:,i) = chance_correlation(pos(i).sep, 100); 
end

figure
for i=1:5
semilogx(C(:,i),'.-')
hold all
%semilogx(M(:,i),'--')
end
axis([0 100 -1 1])

%% isotropy 

for i =1:5
    R(i) = rel_disp(pos(i).sep, ndays, 1);
    iso(:,i) = R(i).avdispzon./R(i).avdispmer;
end

%
figure
plot(iso)
axis([0 100 0 3])

%% velocity correlation
for i=1:5
    Cor(i) = correlation_func(pos(i).sep, ndays);
end

figure
for i = 1:5
    plot((Cor(i).cuu + Cor(i).cvv)/2)
    hold all
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate the Different Metrics %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot_memory_index
plot_chance_correlation
plot_rel_disp_init


%% Dispersion in intermediate regime
% loglog plot 
% Is there a t^3 range? 

figure
loglog(T, dispersion.avdisp/dispersion.avdisp(1),'.-') 
hold all 
loglog(T, 1e-3*T.^3, '--')
%loglog(T, 10^(5.9)*T.^3, '--')

%% semilog plot
% is there an exponential range?

figure
semilogy(T, dispersion.avdisp/dispersion.avdisp(1),'.-') 
hold all 


%% Velocity Correlation

%% Diffusivity 

diff = rel_diff(dispersion.avdisp, ndays); 

% fit a curve to dispersion
fitp = polyfit(T', dispersion.avdisp, 6 ) ;
disp_fit = polyval(fitp, T'); 

% Check fit 
close all
figure
loglog(T, dispersion.avdisp)
hold all
loglog(T, disp_fit)

diff_fit = rel_diff(disp_fit, ndays); 

%% plot diffusivities
Y = dispersion.avdisp.^0.5/1e3;
Yfit = disp_fit.^0.5/1e3;

figure
loglog(Y, diff,'.-') 
hold all 
loglog(Yfit, diff_fit,'o-') 
Ylin = [20:5:100];
loglog(Ylin, Ylin.^2)
loglog(Ylin, 10*Ylin.^(4/3), '--')
axis([10 200 100 10000])

%% Characteristic Time scale tau^2 = <y^2>/2Krel

tau = dispersion.avdisp / 2./ diff; 
tau_fit = disp_fit/2./ diff_fit;

figure
semilogx(Y, tau)
hold all
semilogx(Yfit, tau_fit)
axis([10 300 0 10^7])
