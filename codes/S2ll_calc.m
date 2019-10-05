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


%% obs
[S2ll_obs_deep_mean, S2ll_obs_deep_ci] = S2ll(sep_obs, dist_bin, diff_pres, plevel(2:3)); 
[S2ll_obs_shallow_mean, S2ll_obs_shallow_ci] = S2ll(sep_obs, dist_bin, diff_pres, plevel(1:2)); 

%% model 
S2ll_mod_deep = S2ll_model(mod_sep(2,1).sep, dist_bin); 
S2ll_mod_shallow = S2ll_model(mod_sep(1,1).sep, dist_bin); 

%% 
colors = get(gca, 'ColorOrder'); 
%% Figure for deep
close all 
clear h g 
figure('rend','painters','pos',[10 10 800 600])

h(1) = shadedErrorBar_log(dist_axis/1e3, S2ll_obs_deep_mean, S2ll_obs_deep_ci, ...
            {'-','linewidth',3,'color', colors(1,:)}, 1);
hold all
h(2) = shadedErrorBar_log(dist_axis/1e3, S2ll_obs_shallow_mean, S2ll_obs_shallow_ci, ...
            {'-','linewidth',3,'color', colors(2,:)}, 1);
        
g(1) = loglog(dist_axis/1e3, S2ll_mod_deep, '--', 'linewidth',2,'color', colors(1,:));
g(2) = loglog(dist_axis/1e3, S2ll_mod_shallow, '--', 'linewidth',2,'color', colors(2,:));
%r0 = loglog(dist_axis/1e3, dist_axis./dist_axis*0.1, '-','color', [0.5 0.5 0.5], 'linewidth',1);  
%r23 = loglog(dist_axis/1e3, 100*dist_axis.^(-2/3), '--', 'color', [0.5 0.5 0.5], 'linewidth',1); 
%r2 = loglog(dist_axis/1e3, 1e8*dist_axis.^(-2), '-.', 'color', [0.5 0.5 0.5], 'linewidth',1); 
r23 = loglog(dist_axis/1e3, 1e-6*dist_axis.^(2/3), '-','color', [0.5 0.5 0.5], 'linewidth',1);
r1 =loglog(dist_axis/1e3, 10^(-7.5)*dist_axis.^(1), '--','color', [0.5 0.5 0.5], 'linewidth',1);
r2 =loglog(dist_axis/1e3, 1e-12*dist_axis.^(2), '-.','color', [0.5 0.5 0.5], 'linewidth',1);

A = legend([h(1).mainLine, h(2).mainLine, g, r2, r1, r23], ...
    {'Obs. Deep', 'Obs. Shallow.', 'Mod. Deep', 'Mod. Shallow', '$r^{2}$', '$r^{1}$', '$r^{2/3}$'} ...
    , 'Interpreter','Latex');
legend boxoff
set(A, 'location','best') 
set(gca, 'fontsize', 24) 

axis([1 800 1e-5 2e-2])

xlabel('$r$ (km)', 'Interpreter', 'Latex') 
ylabel('$S2_{ll}$ ($m^{2}/s^{2}$)', 'Interpreter','Latex')

%saveas(gcf,'../figures/S2ll.eps', 'epsc')
%% 
figure 
loglog(dist_axis/1e3, S2ll_obs_deep_mean)
hold all 
loglog(dist_axis/1e3, S2ll_obs_shallow_mean)
loglog(dist_axis/1e3, S2ll_mod_deep)
loglog(dist_axis/1e3, S2ll_mod_shallow)

loglog(dist_axis/1e3, 1e-6*dist_axis.^(2/3))
loglog(dist_axis/1e3, 1e-12*dist_axis.^(2))

axis([1 1000 1e-4 2e-2])


%% 

[S3ll_obs_deep_mean, S3ll_obs_deep_ci] = S3ll(sep_obs, dist_bin, diff_pres, plevel(2:3)); 

%%

figure 
loglog(dist_axis/1e3, abs(S3ll_obs_deep_mean), '.-', 'linewidth',3, 'MarkerSize',20)
hold all 
loglog(dist_axis/1e3, -(S3ll_obs_deep_mean), 'o', 'MarkerSize',15)

loglog(dist_axis/1e3, 10^(-9.5)*dist_axis.^(1), '--','color', [0.5 0.5 0.5], 'linewidth',1);
loglog(dist_axis/1e3, 10^(-19.5)*dist_axis.^(3), '-','color', [0.5 0.5 0.5], 'linewidth',1);

xlabel('$r$ (km)', 'Interpreter', 'Latex') 
ylabel('$S3_{lll}$ ($m^{3}/s^{3}$)', 'Interpreter','Latex')
set(gca,'FontSize',20)

axis([1 1000 1e-7 1e-3])
legend('S3', 'Negative Values', 'r^1', 'r^3', 'Location','Best')



