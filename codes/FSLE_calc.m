% Dhruv Balwada
% 21 September 2018 

% Calculate FSLE 

close all 
clear all 

traj = load ('../data/traj_09_03_13_clean.mat');

% delete sections in the Scotia Sea
traj.Yc(traj.Xc>-82) = NaN;
traj.Xc(traj.Xc>-82) = NaN;

%% Models
d=[4, 10];

for i =1:length(d)
    
    [mod_traj(i).X, mod_traj(i).Y, mod_traj(i).U, mod_traj(i).V, mod_traj(i).T, depth(i)] = loadpairs(d(i));
    
    mod_traj(i).U(mod_traj(i).U==-999) = NaN;
    mod_traj(i).V(mod_traj(i).V==-999) = NaN;
    mod_traj(i).Y(mod_traj(i).X>360-82) = NaN;
    mod_traj(i).U(mod_traj(i).X>360-82) = NaN;
    mod_traj(i).V(mod_traj(i).X>360-82) = NaN;
    mod_traj(i).X(mod_traj(i).X>360-82) = NaN;
    
end

%%
distance_class(1).dist = [5,50]*1e3;
model_sep_calcs

%% interpolate the obs to finer time sections 
Tint = [1:1/24:1117];
Xint = nan(length(Tint),size(traj.Xc,2));
Yint = nan(length(Tint),size(traj.Xc,2));
Uint = nan(length(Tint),size(traj.Xc,2));
Vint = nan(length(Tint),size(traj.Xc,2));
Pint = nan(length(Tint),size(traj.Xc,2));
tempint= nan(length(Tint),size(traj.Xc,2));

for i = 1:size(traj.Xc,2) 
    Xint(:,i) = interp1([1:1117] , traj.Xc(:,i), [1:1/24:1117]);
    Yint(:,i) = interp1([1:1117] , traj.Yc(:,i), [1:1/24:1117]);
    Uint(:,i) = interp1([1:1117] , traj.Uc(:,i), [1:1/24:1117]);
    Vint(:,i) = interp1([1:1117] , traj.Vc(:,i), [1:1/24:1117]);
    Pint(:,i) = interp1([1:1117] , traj.Pi(:,i), [1:1/24:1117]);
    tempint(:,i) = interp1([1:1117] , traj.Ti(:,i), [1:1/24:1117]);
end

traj_int.Xc = Xint; 
traj_int.Yc = Yint;
traj_int.Uc = Uint;
traj_int.Vc = Vint;
traj_int.Pi = Pint;
traj_int.Ti = tempint;
traj_int.name = traj.name;


%% calculate separation time series 

sep_obs_int = calculate_seperation_timeseries(traj_int);
sep_obs = calculate_seperation_timeseries(traj);

%% first define the distance axis 
gamma = 1.4;
dist_bin(1) = 1000; % in m
dist_bin = gamma.^(0:100)*dist_bin(1);
id = find(dist_bin>800*10^3,1);
dist_bin = dist_bin(1:id);
dist_bin(2:end+1) = dist_bin(1:end);
dist_bin(1) = 0;
dist_axis = 0.5*(dist_bin(1:end-1) + dist_bin(2:end));

%%

diff_pres = 100; 
plevel = [500 1000 1800];

flag_1time = 1;
%% 
[fsle_obs_deep_mean, fsle_obs_deep_ci] = fsle(sep_obs, dist_bin, diff_pres, plevel(2:3), 0, flag_1time); 
[fsle_obs_int_deep_mean, fsle_obs_int_deep_ci] = fsle(sep_obs, dist_bin, diff_pres, plevel(2:3), 1, flag_1time); 

[fsle_obs_shallow_mean, fsle_obs_shallow_ci] = fsle(sep_obs, dist_bin, diff_pres, plevel(1:2), 0, flag_1time); 
[fsle_obs_int_shallow_mean, fsle_obs_int_shallow_ci] = fsle(sep_obs, dist_bin, diff_pres, plevel(1:2), 1, flag_1time); 

fsle_mod_deep = fsle_model(mod_sep(2,1).sep, dist_bin, 0, flag_1time); 
fsle_mod_deep_int = fsle_model(mod_sep(2,1).sep, dist_bin, 1, flag_1time); 

fsle_mod_shallow = fsle_model(mod_sep(1,1).sep, dist_bin, 0, flag_1time); 
fsle_mod_shallow_int = fsle_model(mod_sep(1,1).sep, dist_bin, 1, flag_1time); 


%% 
colors = get(gca, 'ColorOrder'); 
%% Figure for deep
close all 
clear h g 
figure('rend','painters','pos',[10 10 800 600])

h(1) = shadedErrorBar_log(dist_axis/1e3, fsle_obs_deep_mean, fsle_obs_deep_ci, ...
            {'--','linewidth',2,'color', colors(1,:)}, 1);
hold all
h(2) = shadedErrorBar_log(dist_axis/1e3, fsle_obs_int_deep_mean, fsle_obs_int_deep_ci, ...
            {'-','linewidth',3,'color', colors(1,:)}, 1);
        
g(1) = loglog(dist_axis/1e3, fsle_mod_deep, '--', 'linewidth',2,'color', colors(2,:));
g(2) = loglog(dist_axis/1e3, fsle_mod_deep_int, '-', 'linewidth',3, 'color', colors(2,:)); 

r0 = loglog(dist_axis/1e3, dist_axis./dist_axis*0.1, '-','color', [0.5 0.5 0.5], 'linewidth',1);  
r23 = loglog(dist_axis/1e3, 100*dist_axis.^(-2/3), '--', 'color', [0.5 0.5 0.5], 'linewidth',1); 
r2 = loglog(dist_axis/1e3, 1e8*dist_axis.^(-2), '-.', 'color', [0.5 0.5 0.5], 'linewidth',1); 


A = legend([h(1).mainLine, h(2).mainLine, g, r0, r23, r2], ...
    {'Obs. Deep', 'Obs. Deep Int.', 'Mod. Deep', 'Mod. Deep Int.', '$r^{0}$', '$r^{-2/3}$', '$r^{-2}$'} ...
    , 'Interpreter','Latex');
legend boxoff
set(A, 'location','best') 
set(gca, 'fontsize', 24) 

axis([1 800 5*1e-3 2])

xlabel('$r$ (km)', 'Interpreter', 'Latex') 
ylabel('$\lambda$ (days $^{-1}$)', 'Interpreter','Latex')

saveas(gcf,'../figures/fsle_deep.eps', 'epsc')

%% Figure for shallow 
close all 
clear h g 
figure('rend','painters','pos',[10 10 800 600])

h(1) = shadedErrorBar_log(dist_axis/1e3, fsle_obs_shallow_mean, fsle_obs_deep_ci, ...
            {'--','linewidth',2,'color', colors(1,:)}, 1);
hold all
h(2) = shadedErrorBar_log(dist_axis/1e3, fsle_obs_int_shallow_mean, fsle_obs_int_deep_ci, ...
            {'-','linewidth',3,'color', colors(1,:)}, 1);
        
g(1) = loglog(dist_axis/1e3, fsle_mod_shallow, '--', 'linewidth',2,'color', colors(2,:));
g(2) = loglog(dist_axis/1e3, fsle_mod_shallow_int, '-', 'linewidth',3, 'color', colors(2,:)); 

r0 = loglog(dist_axis/1e3, dist_axis./dist_axis*0.1, '-','color', [0.5 0.5 0.5], 'linewidth',1);  
r23 = loglog(dist_axis/1e3, 100*dist_axis.^(-2/3), '--', 'color', [0.5 0.5 0.5], 'linewidth',1); 
r2 = loglog(dist_axis/1e3, 1e8*dist_axis.^(-2), '-.', 'color', [0.5 0.5 0.5], 'linewidth',1); 


A = legend([h(1).mainLine, h(2).mainLine, g, r0, r23, r2], ...
    {'Obs. Shallow', 'Obs. Shallow Int.', 'Mod. Shallow', 'Mod. Shallow Int.', '$r^{0}$', '$r^{-2/3}$', '$r^{-2}$'} ...
    , 'Interpreter','Latex');
legend boxoff
set(A, 'location','best') 
set(gca, 'fontsize', 24) 

axis([1 800 5*1e-3 2])

xlabel('$r$ (km)', 'Interpreter', 'Latex') 
ylabel('$\lambda$ (days $^{-1}$)', 'Interpreter','Latex')

saveas(gcf,'../figures/fsle_shallow.eps', 'epsc')

%%  Calculate the zonal and meridional FSLE seperately. 
% turns out that the system is more or less isotropic in this. 

[fsle_obs_zon, fsle_obs_mer] = fsle_zonmer(sep_obs, dist_bin, diff_pres, plevel(2:3), 0, flag_1time); 

%% 
figure 
loglog(dist_axis/1e3 , fsle_obs_zon) 
hold all 
loglog(dist_axis/1e3 , fsle_obs_mer) 