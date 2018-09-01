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

%% Number of pairs as a function time 

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
legend('500-1000m', '1000-1800m')
xlabel('Days')
ylabel('No. Pairs')
set(gca,'fontsize',20)

saveas(gca,'../figures/no_pairs.pdf')


%% Initial distribution of pairs
for i = 1:pairs_ncon 
   ini_dist(i) = pos_ncon.sep(i).dist(1); 
end
    
%
figure
hist(ini_dist)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate the Different Metrics %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Memory Index (position Correlation)
% Equation 5.1 in Foussard et al 2017 
% M(t) = <Y.Y_0>/ Y_0<Y^2>^0.5
% Where Y is the separation vector, and Y_o is the initial separation
% vector

for i = 1:500
    y = randsample(pairs_ncon,pairs_ncon,true);
    
    M(1:ndays,i) = memory_index(pos_ncon.sep(y), 100);
end

%%
close all
figure
semilogx(M, 'linewidth', 0.5)
hold all

semilogx(nanmean(M,2), 'linewidth',2)

grid

%% Correlation to initial position <dv.y_o> (should be zero for original pairs)

for i = 1:500
    y = randsample(pairs_ncon,pairs_ncon,true);

    C(1:ndays,i) = chance_correlation(pos_ncon.sep(y), 100);
end

%C = chance_correlation(pos_ncon.sep, 100);

%%
figure
plot(C)
grid

%% Dispersion during ballistic regime 
for i = 1:500
    y = randsample(pairs_ncon,pairs_ncon,true);

    [dispersion(i)] = rel_disp(pos_ncon.sep,ndays,23);
end

%% 
figure
loglog(T, dispersion(1).avdisp-dispersion(1).avdisp(1),'.-') 
hold all 
loglog(T, 10^(7.3)*T.^2, '--')
%loglog(T, 10^(5.9)*T.^3, '--')

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

%% Isotropy

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
