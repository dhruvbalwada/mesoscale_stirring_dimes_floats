% Dhruv Balwada
% Created 8 Jan 2020

% influence of NIWs on metrics that bin in space
% 
% one of the problems of this test below is that we had to choose A to be
% random number. One could envision a more complicated model in which there
% is some spatial coherence 

clear all 
close all
% load some model data 
d = 10; 

distance = [10,15]*1e3;

[mod_traj.X, mod_traj.Y, mod_traj.U, mod_traj.V, mod_traj.T, depth] = loadpairs(d);
    
mod_traj.U(mod_traj.U==-999) = NaN;
mod_traj.V(mod_traj.V==-999) = NaN;
mod_traj.Y(mod_traj.X>360-70) = NaN;
mod_traj.U(mod_traj.X>360-70) = NaN;
mod_traj.V(mod_traj.X>360-70) = NaN;
mod_traj.X(mod_traj.X>360-70) = NaN;

%%
f = 2*2*pi/(24*3600)*cosd(-55);
%%
mod_sep_orig = model_sep_calcs_w_NIWs(mod_traj, distance, f, 0,0);
mod_sep_1 = model_sep_calcs_w_NIWs(mod_traj, distance, f, 0,0.3);
mod_sep_2 = model_sep_calcs_w_NIWs(mod_traj, distance, f, 0,0.6);
mod_sep_3 = model_sep_calcs_w_NIWs(mod_traj, distance, f, 0,1);
mod_sep_4 = model_sep_calcs_w_NIWs(mod_traj, distance, f, 0,2);


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
fsle_mod_1 = fsle_model(mod_sep_1.sep, dist_bin, 1, flag_1time); 
fsle_mod_2 = fsle_model(mod_sep_2.sep, dist_bin, 1, flag_1time); 
fsle_mod_3 = fsle_model(mod_sep_3.sep, dist_bin, 1, flag_1time); 
fsle_mod_4 = fsle_model(mod_sep_4.sep, dist_bin, 1, flag_1time); 
%fsle_mod_5 = fsle_model(mod_sep_5.sep, dist_bin, 1, flag_1time); 
%fsle_mod_6 = fsle_model(mod_sep_6.sep, dist_bin, 1, flag_1time); 


%%
figure
loglog(dist_axis/1e3, fsle_mod_orig, '-.', 'linewidth',2);
hold all 
%loglog(dist_axis/1e3, fsle_mod_1, '-', 'linewidth',2);
%loglog(dist_axis/1e3, fsle_mod_2, '-', 'linewidth',2);
loglog(dist_axis/1e3, fsle_mod_3, '-', 'linewidth',2);
%loglog(dist_axis/1e3, fsle_mod_4, '--', 'linewidth',2);
%loglog(dist_axis/1e3, fsle_mod_5, '--', 'linewidth',2);
%loglog(dist_axis/1e3, fsle_mod_6, '--', 'linewidth',2);

loglog(dist_axis/1e3, 100*dist_axis.^(-2/3))
loglog(dist_axis/1e3, 1000*dist_axis.^(-1))

%% S2ll
