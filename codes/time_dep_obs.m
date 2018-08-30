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
traj.Yc(traj.Xc>-70) = NaN;
traj.Xc(traj.Xc>-70) = NaN;

% calculate the separation time series for all possible pairs.
sep = calculate_seperation_timeseries(traj);

%% Specify pressure criterion, minimum days, and distance classes.
plevel = [1000 1800];
ndays = 50;
pdiff =150;

%% 
T = 0:ndays-1;

%% Find pairs corresponding to the criterion.
distance_class_obs = [20,25]*1000;

[pairs_ncon,pos_ncon] = find_pairs(sep, distance_class_obs, plevel, pdiff, ndays);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate the Different Metrics %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Memory Index (position Correlation)
% Equation 5.1 in Foussard et al 2017 
% M(t) = <Y.Y_0>/ Y_0<Y^2>^0.5
% Where Y is the separation vector, and Y_o is the initial separation
% vector

M = memory_index(pos_ncon.sep, 100);

%%
close all
figure
semilogx(M)
grid

%% Correlation to initial position <dv.y_o> (should be zero for original pairs)

C = chance_correlation(pos_ncon.sep, 100);

%%
figure
plot(C)
grid

%% Dispersion during ballistic regime 

[dispersion] = rel_disp(pos_ncon.sep,ndays,23);

%% 
figure
loglog(T, dispersion.avdisp-dispersion.avdisp(1),'.-') 
hold all 
loglog(T, 10^(7.3)*T.^2, '--')
loglog(T, 10^(5.9)*T.^3, '--')

%% Dispersion in intermediate regime
%% loglog plot 
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


