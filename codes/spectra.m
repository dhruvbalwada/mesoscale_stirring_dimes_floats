%% Calculate the rotary spectra for the model and RAFOS floats

clear all
close all 

%% load data 

obs = load ('../data/traj_cut_DIMES_120days.mat');

%%
Xm = nanmean(obs.X_pday,1);
Pm = nanmean(obs.Pi_pday,1);
id = find(Pm>500 & Pm<1000 & Xm<-80); 
%%
d = 3;
for i =1:length(d)
    
    [mod_traj(i).X, mod_traj(i).Y, mod_traj(i).U, mod_traj(i).V, mod_traj(i).T, depth(i)] = loadpairs(d(i));
    
    mod_traj(i).U(mod_traj(i).U==-999) = NaN;
    mod_traj(i).V(mod_traj(i).V==-999) = NaN;
    mod_traj(i).Y(mod_traj(i).X>360-70) = NaN;
    mod_traj(i).U(mod_traj(i).X>360-70) = NaN;
    mod_traj(i).V(mod_traj(i).X>360-70) = NaN;
    mod_traj(i).X(mod_traj(i).X>360-70) = NaN;
    
end

%% Rotary spec using Jlab
X =  obs.X_pday;
nsegs = size(X,2); 
len = size(X,1); 
psi = sleptap(len); 


CV_obs = zeros(len,nsegs);
% complex velocities
for nseg = 1:nsegs
    CV_obs(:,nseg) = latlon2uv(obs.T_pday(:,nseg), obs.Y_pday(:,nseg), obs.X_pday(:,nseg))/100;
    CV_obs_forward(:,nseg) = latlon2uv(obs.T_pday(:,nseg), obs.Y_pday(:,nseg), obs.X_pday(:,nseg), 'forward')/100;
end

CV_mod = zeros(len, 1200); 
for nseg=1:1200
    CV_mod(:,nseg) = mod_traj.U(2:121, nseg) + sqrt(-1)*mod_traj.V(2:121, nseg);
    CV_mod2(:,nseg) = latlon2uv(mod_traj.T(2:121), mod_traj.Y(2:121,nseg), mod_traj.X(2:121,nseg) )/100;
    CV_mod3(:,nseg) = latlon2uv(mod_traj.T(2:121), mod_traj.Y(2:121,nseg), mod_traj.X(2:121,nseg), 'forward' )/100;
end
    
%%
% Spectrum

[F_obs, SPP_obs, SNN_obs, SPN_obs] = mspec(CV_obs, psi,'cyclic');
[F_obs_forward, SPP_obs_forward, SNN_obs_forward, SPN_obs_forward] = mspec(CV_obs_forward, psi,'cyclic');
[F_mod, SPP_mod, SNN_mod, SPN_mod] = mspec(CV_mod, psi,'cyclic');
[F_mod2, SPP_mod2, SNN_mod2, SPN_mod2] = mspec(CV_mod2, psi,'cyclic');
[F_mod3, SPP_mod3, SNN_mod3, SPN_mod3] = mspec(CV_mod3, psi,'cyclic');

%F = F;


%% 
figure 
loglog(F_obs, mean(SPP_obs_forward(:, id),2),'--', 'linewidth',2,'color', [0.3010    0.7450    0.9330])
%loglog(F_obs, mean(SPP_obs(:, id),2),'--', 'linewidth',2,'color', [0.3010    0.7450    0.9330])
hold all
%loglog(F_obs, mean(SNN_obs(:, id),2), '-', 'linewidth',2,'color', [0.3010    0.7450    0.9330])
loglog(F_obs, mean(SNN_obs_forward(:, id),2), '-', 'linewidth',2,'color', [0.3010    0.7450    0.9330])

%loglog(F_mod, nanmean(SPP_mod,2), 'linewidth',1, 'color',[0.8500    0.3250    0.0980])
%loglog(F_mod, nanmean(SNN_mod,2), '--', 'linewidth',1, 'color',[0.8500    0.3250    0.0980])

%loglog(F_mod, nanmean(SPP_mod2,2), '--','linewidth',2, 'color',[0.8500    0.3250    0.0980])
loglog(F_mod, nanmean(SPP_mod3,2), '--','linewidth',2, 'color',[0.8500    0.3250    0.0980])
loglog(F_mod, nanmean(SNN_mod3,2), '-', 'linewidth',2, 'color',[0.8500    0.3250    0.0980])

loglog(F_obs, 3e-6*F_obs.^-3, '-', 'color',[0.5 0.5 0.5],'linewidth',1.5)
loglog(F_obs, 1e-9*F_obs.^-5, '--', 'color',[0.5 0.5 0.5],'linewidth',1.5)
%loglog(f*[1,1] , [0.1,100])
xlim([1/120, 1/2])
%ylim([5e-6, 0.1])
ylim([5e-6, 0.5])

legend('Obs. Anticyclonic', 'Obs. Cyclonic', 'Mod. Anticyclonic', ...
    'Mod. Cyclonic','\omega^{-3}', '\omega^{-5}','location','best')
set(gca,'fontsize',20)
xlabel('\omega (cycles/day)')
ylabel('KE Spectral Density (m^2/s)')
    
    