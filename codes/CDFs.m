% CDFs model and obs 

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

%%
mod_sep = model_sep_calcs(mod_traj, distance); 
%%
disp = rel_disp(mod_sep.sep, 100);

%% 
figure
plot(disp.disp.^0.5)
%% 
figure
%histogram(disp.disp(1,:).^0.5/1e3, 'Normalization','pdf')
hold all 
histogram(disp.disp(5,:).^0.5/1e3, 'Normalization','pdf')
histogram(disp.disp(10,:).^0.5/1e3, 'Normalization','pdf')
histogram(disp.disp(15,:).^0.5/1e3, 'Normalization','pdf')

%%
edges = linspace(0, 250, 170); 
%% 
figure
%histogram(disp.disp(2,:).^0.5/1e3, 'Normalization','pdf')
hold all 
h(1) = histogram(disp.disp(5,:).^0.5/1e3, edges, 'Normalization','cdf'); 
h(2) = histogram(disp.disp(10,:).^0.5/1e3, edges,  'Normalization','cdf');
h(3) = histogram(disp.disp(15,:).^0.5/1e3, edges,  'Normalization','cdf');
h(4) = histogram(disp.disp(20,:).^0.5/1e3, edges,  'Normalization','cdf');
h(5) = histogram(disp.disp(25,:).^0.5/1e3, edges,  'Normalization','cdf');
grid on

%%
clear sep_cdf
for i =1:100
    sep_cdf(:,i) =  histcounts(disp.disp(i,:).^0.5/1e3, edges, 'Normalization','cdf'); 
end

%%
T = 1:100;    
r = 0.5*(edges(1:end-1) + edges(2:end));
%% 

figure 
contourf(r, T, sep_cdf', 'Edgecolor','none')
hold all
contour(r, T, sep_cdf', [0.1, 0.9], 'Edgecolor', 'k','linewidth',2)
plot(disp.avdisp.^0.5/1e3, T, 'Linewidth',2) 

    