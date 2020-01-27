% Created 17 September 2018
% Dhruv Balwada

clear all
close all

%% Load data
% Observations
obs_traj = load ('../data/traj_09_03_13_clean.mat');

% delete sections in the Scotia Sea.
obs_traj.Yc(obs_traj.Xc>-82) = NaN;
obs_traj.Xc(obs_traj.Xc>-82) = NaN;

% Models
d=[4, 10];

for i =1:length(d)
    [mod_traj(i).X, mod_traj(i).Y, mod_traj(i).U, mod_traj(i).V, mod_traj(i).T, depth(i)] = loadpairs(d(i));
    
    mod_traj(i).U(mod_traj(i).U==-999) = NaN;
    mod_traj(i).V(mod_traj(i).V==-999) = NaN;
    mod_traj(i).Y(mod_traj(i).X>360-70) = NaN;
    mod_traj(i).U(mod_traj(i).X>360-70) = NaN;
    mod_traj(i).V(mod_traj(i).X>360-70) = NaN;
    mod_traj(i).X(mod_traj(i).X>360-70) = NaN;
    
end

%%
% calculate the separation time series for all possible pairs.
% Observations are different from model in how pairs are defined.
% The separation in pressure is considered at a later stage for obs, but is
% initially considered in model.

sep_obs = calculate_seperation_timeseries(obs_traj);

%% consider 3 distance classes 10-15, 30-35, 50-55.
clear distance_class obs_traj
distance_class(1).dist = [10,15]*1e3;
distance_class(2).dist = [20,25]*1e3;
distance_class(3).dist = [30,35]*1e3;
distance_class(4).dist = [40,45]*1e3;
distance_class(5).dist = [50,55]*1e3;
distance_class(6).dist = [60,65]*1e3;

%%
plevel = [500 1000 1800];
ndays  = 100;
pdiff  = 200;

%%
T=[0:ndays-1];

%%
for i =1:length(distance_class)
    for j=1:length(plevel)-1
        [obs_pairs(j,i), obs_sep(j,i)] = find_pairs(sep_obs, distance_class(i).dist, [plevel(j) plevel(j+1)], pdiff, ndays);
    end
end

%% calc sep models
for j = 1:length(distance_class)
    for i =1:length(mod_traj)
        mod_sep(j,i) = model_sep_calcs(mod_traj(i), distance_class(j).dist); 
    end
end
%%
for i =1:length(distance_class)
    for j=1:2
        mod_pairs(j,i) = length(mod_sep(j,i).sep);
    end
end

clear mod_traj
%%
mod_dist_ini = zeros(2,length(distance_class)); 
dist = nan*ones(1176, 2,length(distance_class));

for i =1:length(distance_class)
    for j=1:2
        for k =1:length(mod_sep(j,i).sep)
            mod_dist_ini(j,i) = mod_dist_ini(j,i) + mod_sep(j,i).sep(k).dist(1); 
            dist(k,j,i) = mod_sep(j,i).sep(k).dist(1);
        end
    end
end
mod_dist_ini = mod_dist_ini./mod_pairs; 

%% Things that can be calc after this : 
% plot_memory_index_all.m
% plot_chance_correlation_all.m 
