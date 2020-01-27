% Dhruv Balwada
% 21 September 2018 

% Calculate S2ll 

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

%% calculate separation time series 

sep_obs = calculate_seperation_timeseries(traj);

%% first define the distance axis 
gamma = 1.45;
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


%% obs

rd_obs_deep_struct = rel_diff_inst(sep_obs, dist_bin, diff_pres, plevel(2:3),2); 
rd_obs_shallow_struct = rel_diff_inst(sep_obs, dist_bin, diff_pres, plevel(1:2),2); 

%[S2ll_obs_deep_mean, S2ll_obs_deep_ci] = S2ll(sep_obs, dist_bin, diff_pres, plevel(2:3)); 
%[S2ll_obs_shallow_mean, S2ll_obs_shallow_ci] = S2ll(sep_obs, dist_bin, diff_pres, plevel(1:2)); 

%% model 
rd_mod_deep = rel_diff_inst_model(mod_sep(2,1).sep, dist_bin,2); 
rd_mod_shallow = rel_diff_inst_model(mod_sep(1,1).sep, dist_bin,2); 

%%
colors = get(gca, 'ColorOrder'); 

%% deep only 

close all 
clear h g 
figure('rend','painters','pos',[10 10 800 600])

%h(1) = shadedErrorBar_log(dist_axis/1e3, rd_obs_deep_mean, rd_obs_deep_ci, ...
%            {'-','linewidth',3,'color', colors(1,:)}, 1);
h(1) = shadedErrorBar_log(dist_axis/1e3, rd_obs_deep_struct.rdmean, rd_obs_deep_struct.rdci, ...
            {'-','linewidth',3,'color', colors(1,:)}, 1);

hold all 
h(2) = shadedErrorBar_log(dist_axis/1e3, rd_obs_deep_struct.rdXmean, rd_obs_deep_struct.rdXci, ...
            {'--','linewidth',2,'color', colors(1,:)}, 1);

h(3) = shadedErrorBar_log(dist_axis/1e3, rd_obs_deep_struct.rdYmean, rd_obs_deep_struct.rdYci, ...
            {'-.','linewidth',2,'color', colors(1,:)}, 1);

r2 = loglog(dist_axis/1e3, 1e-6*dist_axis.^2, '-.','color', [0.5 0.5 0.5], 'linewidth',1);
r43 = loglog(dist_axis/1e3, 4e-4*dist_axis.^(4/3), '--','color', [0.5 0.5 0.5], 'linewidth',1);

x= [500, 500]; 
y = 2*[800,1200];
kbal = loglog(x, y, 'o-', 'linewidth',2);

A = legend([h(1).mainLine, h(2).mainLine, h(3).mainLine,r2, r43, kbal], ...
    {'K','$K_x$','$K_y$','$r^2$','$r^{4/3}$','$2K_{single}$'} ...
    , 'Interpreter','Latex');
legend boxoff
set(A, 'location','best') 
set(gca, 'fontsize', 24) 

axis([1 800 1 2e4])

xlabel('$r$ (km)', 'Interpreter', 'Latex') 
ylabel('$K (r)$ ($m^{2}/s$)', 'Interpreter','Latex')

%% shallow only 

%close all 
clear h g 
figure('rend','painters','pos',[10 10 800 600])

%h(1) = shadedErrorBar_log(dist_axis/1e3, rd_obs_deep_mean, rd_obs_deep_ci, ...
%            {'-','linewidth',3,'color', colors(1,:)}, 1);
h(1) = shadedErrorBar_log(dist_axis/1e3, rd_obs_shallow_struct.rdmean, rd_obs_shallow_struct.rdci, ...
            {'-','linewidth',3,'color', colors(1,:)}, 1);

hold all 
h(2) = shadedErrorBar_log(dist_axis/1e3, rd_obs_shallow_struct.rdXmean, rd_obs_shallow_struct.rdXci, ...
            {'--','linewidth',2,'color', colors(1,:)}, 1);

h(3) = shadedErrorBar_log(dist_axis/1e3, rd_obs_shallow_struct.rdYmean, rd_obs_shallow_struct.rdYci, ...
            {'-.','linewidth',2,'color', colors(1,:)}, 1);

r2 = loglog(dist_axis/1e3, 3e-6*dist_axis.^2, '-.','color', [0.5 0.5 0.5], 'linewidth',1);
r43 = loglog(dist_axis/1e3, 9e-4*dist_axis.^(4/3), '--','color', [0.5 0.5 0.5], 'linewidth',1);

x= [500, 500]; 
y = 2*[550,850];
kbal = loglog(x, y, 'o-','linewidth',2);

A = legend([h(1).mainLine, h(2).mainLine, h(3).mainLine,r2, r43, kbal], ...
    {'K','$K_x$','$K_y$','$r^2$','$r^{4/3}$','$2K_{single}$'} ...
    , 'Interpreter','Latex');
legend boxoff
set(A, 'location','best') 
set(gca, 'fontsize', 24) 

axis([1 800 1 2e4])

xlabel('$r$ (km)', 'Interpreter', 'Latex') 
ylabel('$K (r)$ ($m^{2}/s$)', 'Interpreter','Latex')


%% model deep 

clear h g 
figure('rend','painters','pos',[10 10 800 600])
g(1) = loglog(dist_axis/1e3, rd_mod_deep.rdmean,'-', 'linewidth',3,'color', colors(2,:));
hold all
g(2) = loglog(dist_axis/1e3, rd_mod_deep.rdXmean,'--', 'linewidth',2,'color', colors(2,:));
g(3) = loglog(dist_axis/1e3, rd_mod_deep.rdYmean,'-.', 'linewidth',2,'color', colors(2,:));

r2 = loglog(dist_axis/1e3, 1e-6*dist_axis.^2, '-.','color', [0.5 0.5 0.5], 'linewidth',1);
r43 = loglog(dist_axis/1e3, 9e-4*dist_axis.^(4/3), '--','color', [0.5 0.5 0.5], 'linewidth',1);

x= [500, 500]; 
y = 2*[800,1200];
kbal = loglog(x, y, 'o-', 'linewidth',2);

A = legend([g(1), g(2), g(3),r2, r43, kbal], ...
    {'K','$K_x$','$K_y$','$r^2$','$r^{4/3}$','$2K_{single}$'} ...
    , 'Interpreter','Latex');
legend boxoff
set(A, 'location','best') 
set(gca, 'fontsize', 24) 

axis([1 800 1 2e4])

xlabel('$r$ (km)', 'Interpreter', 'Latex') 
ylabel('$K (r)$ ($m^{2}/s$)', 'Interpreter','Latex')

%% model shallow 

clear h g 
figure('rend','painters','pos',[10 10 800 600])
g(1) = loglog(dist_axis/1e3, rd_mod_shallow.rdmean,'-', 'linewidth',3,'color', colors(2,:));
hold all
g(2) = loglog(dist_axis/1e3, rd_mod_shallow.rdXmean,'--', 'linewidth',2,'color', colors(2,:));
g(3) = loglog(dist_axis/1e3, rd_mod_shallow.rdYmean,'-.', 'linewidth',2,'color', colors(2,:));

r2 = loglog(dist_axis/1e3, 1e-6*dist_axis.^2, '-.','color', [0.5 0.5 0.5], 'linewidth',1);
r43 = loglog(dist_axis/1e3, 9e-4*dist_axis.^(4/3), '--','color', [0.5 0.5 0.5], 'linewidth',1);

x= [500, 500]; 
y = 2*[550,850];
kbal = loglog(x, y, 'o-', 'linewidth',2);

A = legend([g(1), g(2), g(3),r2, r43, kbal], ...
    {'K','$K_x$','$K_y$','$r^2$','$r^{4/3}$','$2K_{single}$'} ...
    , 'Interpreter','Latex');
legend boxoff
set(A, 'location','best') 
set(gca, 'fontsize', 24) 

axis([1 800 1 2e4])

xlabel('$r$ (km)', 'Interpreter', 'Latex') 
ylabel('$K (r)$ ($m^{2}/s$)', 'Interpreter','Latex')



