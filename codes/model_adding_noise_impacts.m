% Dhruv Balwada
% Created 8 Jan 2020

% influence of noise on metrics that bin in space
% noise - internal waves or instrument tracking noise. 
% 
clear all 
close all
% load some model data 
d = 10; 

distance = [10,15]*1e3;

i = 1; 
[mod_traj(i).X, mod_traj(i).Y, mod_traj(i).U, mod_traj(i).V, mod_traj(i).T, depth(i)] = loadpairs(d(i));
    
mod_traj(i).U(mod_traj(i).U==-999) = NaN;
mod_traj(i).V(mod_traj(i).V==-999) = NaN;
mod_traj(i).Y(mod_traj(i).X>360-70) = NaN;
mod_traj(i).U(mod_traj(i).X>360-70) = NaN;
mod_traj(i).V(mod_traj(i).X>360-70) = NaN;
mod_traj(i).X(mod_traj(i).X>360-70) = NaN;

%% add some noise 

Lk = 3/111; % noise in [km] needs to be scaled to represent lon or lat

xrand = randn(736, 1200);
yrand = randn(736, 1200);

mod_traj_orig = mod_traj;

mod_traj.X = mod_traj.X + xrand*Lk.*cosd(mod_traj.Y);
mod_traj.Y = mod_traj.Y + yrand*Lk;

%%
figure 
hold all 
for i = 1:3
    plot(mod_traj_orig.X(:,i), mod_traj_orig.Y(:,i))
    plot(mod_traj.X(:,i), mod_traj.Y(:,i), '--')
end

%% 
mod_sep_orig = model_sep_calcs(mod_traj_orig, distance); 
mod_sep_rand = model_sep_calcs(mod_traj, distance); 

%% first define the distance axis 
gamma = 1.4;
dist_bin(1) = 1000; % in m
dist_bin = gamma.^(0:100)*dist_bin(1);
id = find(dist_bin>800*10^3,1);
dist_bin = dist_bin(1:id);
dist_bin(2:end+1) = dist_bin(1:end);
dist_bin(1) = 0;
dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));

flag_1time = 1;
%% FSLE
fsle_mod_orig = fsle_model(mod_sep_orig.sep, dist_bin, 1, flag_1time); 
fsle_mod_rand = fsle_model(mod_sep_rand.sep, dist_bin, 1, flag_1time); 

%%
figure
loglog(dist_axis/1e3, fsle_mod_orig, '--', 'linewidth',2);
hold all
loglog(dist_axis/1e3, fsle_mod_rand, '--', 'linewidth',2);

%% relative diff

rd_mod_orig = rel_diff_inst_model(mod_sep_orig.sep, dist_bin, 3); 
rd_mod_rand = rel_diff_inst_model(mod_sep_rand.sep, dist_bin, 3); 

%%
figure
loglog(dist_axis/1e3, rd_mod_orig.rdmean, '--', 'linewidth',2);
hold all
loglog(dist_axis/1e3, rd_mod_rand.rdmean, '--', 'linewidth',2);