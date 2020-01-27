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
X =mod_traj.X; 
Y =mod_traj.Y; 
%% Interpolate to hourly
T = 0:735; 
Tint = 0:1/24:735;
   
Xint = interp1(T, X, Tint, 'linear'); 
Yint = interp1(T, Y, Tint, 'linear'); 

%%
figure
plot(T, X)
hold all
plot(Tint, Xint)

%% adding NIWs
A = 3*randn(1200, 1);
A = abs(A);

f = 2* (2*pi/24/3600) * sind(-55); 
%%
Xniw = 0*Xint; 
Yniw = 0*Yint; 

for i = 1:1200
    Xw = A(i)*sin(f*Tint*24*3600);
    Yw = A(i)*(cos(f*Tint*24*3600)-1);
    
    [Ylw, Xlw] = xy2latlon(Xw, Yw, -55, 270);
    
    Xniw(:,i) = Xint(:,i) + Xlw' - Xlw(1); 
    Yniw(:,i) = Yint(:,i) + Ylw' - Ylw(1); 
end

%% 
i =1:3; 
l = 150; 
figure, 
plot(Xint(1:l, i), Yint(1:l, i)) 
hold all
plot(Xint(1, i), Yint(1, i), '*') 
plot(Xniw(1:l, i), Yniw(1:l, i)) 
plot(Xniw(1, i), Yniw(1, i), 'o') 

%% Resample
Xrs = interp1(Tint, Xniw, T, 'nearest');
Yrs = interp1(Tint, Yniw, T, 'nearest');

%% 
len = 120; 
CV_mod = zeros(len, 1200); 
CV_mod2 = CV_mod;
%CV_mod3 = CV_mod;
clear CV_mod3

for nseg=1:1200
    CV_mod(:,nseg) = latlon2uv(mod_traj.T(2:121), mod_traj.Y(2:121,nseg), mod_traj.X(2:121,nseg) )/100;
    CV_mod2(:,nseg) = latlon2uv(mod_traj.T(2:121), Yrs(2:121,nseg), Xrs(2:121,nseg) )/100;
    CV_mod3(:,nseg) = latlon2uv(Tint(2:121*24), Yint(2:121*24,nseg), Xint(2:121*24,nseg) )/100;
    CV_mod4(:,nseg) = latlon2uv(Tint(2:121*24), Yniw(2:121*24,nseg), Xniw(2:121*24,nseg) )/100;
end

%% 

psi = sleptap(len); 
psi2 = sleptap(2903);

%%
[F_mod, SPP_mod, SNN_mod, SPN_mod] = mspec(CV_mod, psi,'cyclic');
[F_mod2, SPP_mod2, SNN_mod2, SPN_mod2] = mspec(CV_mod2, psi,'cyclic');
%%
[F_mod3, SPP_mod3, SNN_mod3, SPN_mod3] = mspec(CV_mod3, psi2,'cyclic');
[F_mod4, SPP_mod4, SNN_mod4, SPN_mod4] = mspec(CV_mod4, psi2,'cyclic');

%%
figure
loglog(F_mod3, nanmean(SPP_mod3,2))
hold all 
loglog(F_mod3, nanmean(SNN_mod3,2))

loglog(F_mod3, nanmean(SPP_mod4,2))
loglog(F_mod3, nanmean(SNN_mod4,2))

%%
figure 
loglog(F_mod3, SPP_mod4(:,1))
hold all 
loglog(F_mod3, SNN_mod4(:,1))
%% 
figure 

loglog(F_mod, nanmean(SPP_mod,2), 'linewidth',1, 'color',[0.3010    0.7450    0.9330])
hold all
loglog(F_mod, nanmean(SNN_mod,2), '--', 'linewidth',1, 'color',[0.3010    0.7450    0.9330])

loglog(F_mod, nanmean(SPP_mod2,2), '--','linewidth',2, 'color',[0.8500    0.3250    0.0980])
loglog(F_mod, nanmean(SNN_mod2,2), '-', 'linewidth',2, 'color',[0.8500    0.3250    0.0980])

%loglog(F_obs, 1e-5*F_obs.^-3, '-', 'color',[0.5 0.5 0.5],'linewidth',1.5)
%loglog(F_obs, 1e-8*F_obs.^-5, '--', 'color',[0.5 0.5 0.5],'linewidth',1.5)
%loglog(f*[1,1] , [0.1,100])
xlim([1/240, 1])
ylim([5e-6, 0.3])
legend('Obs. Anticyclonic', 'Obs. Cyclonic', 'Mod. Anticyclonic', ...
    'Mod. Cyclonic','\omega^{-3}', '\omega^{-5}','location','best')
set(gca,'fontsize',20)
xlabel('\omega (cycles/day)')
ylabel('KE Spectral Density (m^2/s)')

%% 

mod_traj_orig = mod_traj; 
mod_traj_niw = mod_traj; 
mod_traj_niw.X = Xrs;
mod_traj_niw.Y = Yrs;

%% 
CV_rs = latlon2uv(T, Yrs, Xrs);
mod_traj_niw.U = real(CV_rs);
mod_traj_niw.V = imag(CV_rs);

%% 
mod_sep_orig = model_sep_calcs(mod_traj_orig, distance); 
mod_sep_niw = model_sep_calcs(mod_traj_niw, distance); 

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
fsle_mod_niw = fsle_model(mod_sep_niw.sep, dist_bin, 1, flag_1time); 

%%
figure
loglog(dist_axis/1e3, fsle_mod_orig, '--', 'linewidth',2);
hold all
loglog(dist_axis/1e3, fsle_mod_niw, '--', 'linewidth',2);


%% relative diff

rd_mod_orig = rel_diff_inst_model(mod_sep_orig.sep, dist_bin, 2); 
rd_mod_niw = rel_diff_inst_model(mod_sep_niw.sep, dist_bin, 2); 

%%
figure
loglog(dist_axis/1e3, rd_mod_orig.rdmean, '--', 'linewidth',2);
hold all
loglog(dist_axis/1e3, rd_mod_niw.rdmean, '--', 'linewidth',2);
loglog(dist_axis/1e3, 1e-7*dist_axis.^2)
loglog(dist_axis/1e3, 1e-7*dist_axis.^(4/3))