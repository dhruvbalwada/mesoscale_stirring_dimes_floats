% Master function that does all the fitting for the PDFs
% change these for different plots

% run main_obs_model.m
% then run plot_rel_disp_all.m

%% supplementary functions
flung = @(Tlung,xdata, Ro)(Ro^2)*exp(8*xdata/Tlung);
frich = @(bguess,xdata, Ro) factorial(5)/2* ...
    ((4*bguess*xdata/9).^3).* ...
    exp(-9*Ro^(2/3)/4/bguess./xdata).* ...
    hypergeom(6,3,9*Ro^(2/3)/4/bguess./xdata);
fdiff = @(kappa, xdata, Ro) Ro^2 + 4*kappa*xdata;

%%
clear fit_parms_mod fit_parms_obs
a=5;
for j=1:2
    dispersion = Obs_disp(j,1).avdisp;
    fit_parms_obs(j) = dispersion_fitting(dispersion, T', a);
    fit_parms_obs_log(j) = dispersion_fitting_log(dispersion, T', a);
end

for j=1:2
    dispersion = Mod_disp(j,1).avdisp;
    fit_parms_mod(j) = dispersion_fitting(dispersion, T', a);
    fit_parms_mod_log(j) = dispersion_fitting_log(dispersion, T', a);
end

%%
colors = get(gca, 'ColorOrder'); 

%% Plots of dispersions and fits 

% Obs 

j=2; % shallow (1), deep (2)
figure('rend','painters','pos',[10 10 400 300])
loglog(T, Obs_disp(j,1).avdisp, 'linewidth',3, 'color', colors(4-j,:))
hold all
loglog(T, flung(fit_parms_obs(j).Tlung, T, fit_parms_obs(j).Ro), 'color',[0.5 0.5 0.5], 'linewidth',1)
loglog(T, frich(fit_parms_obs(j).beta , T, fit_parms_obs(j).Ro),'--','color',[0.5 0.5 0.5], 'linewidth',1)
loglog(T, fdiff(fit_parms_obs(j).kappa , T, fit_parms_obs(j).Ro),'-.','color',[0.5 0.5 0.5], 'linewidth',1)


% loglog(T, flung(fit_parms_obs_log(j).Tlung, T, fit_parms_obs_log(j).Ro), 'color','k', 'linewidth',1)
% loglog(T, frich(fit_parms_obs_log(j).beta , T, fit_parms_obs_log(j).Ro),'--','color','k', 'linewidth',1)
% loglog(T, fdiff(fit_parms_obs_log(j).kappa , T, fit_parms_obs_log(j).Ro),'-.','color','k', 'linewidth',1)

axis([1 30 1e8 1e10])
xlabel('t (Days)')
ylabel('D^2 (t) (m^2)')
set(gca, 'fontsize',18)


%if j==1
%    saveas(gcf,'../figures/fit_disp_obs_shallow.eps', 'epsc')
%else
%    saveas(gcf,'../figures/fit_disp_obs_deep.eps', 'epsc')
%end

%%
% Mod 
j=2; 
figure('rend','painters','pos',[10 10 400 300])
loglog(T, Mod_disp(j,1).avdisp, 'linewidth',3, 'color', colors(4-j,:))
hold all
loglog(T, flung(fit_parms_mod(j).Tlung, T, fit_parms_mod(j).Ro), 'color',[0.5 0.5 0.5], 'linewidth',1)
loglog(T, frich(fit_parms_mod(j).beta , T, fit_parms_mod(j).Ro),'--','color',[0.5 0.5 0.5], 'linewidth',1)
loglog(T, fdiff(fit_parms_mod(j).kappa , T, fit_parms_mod(j).Ro),'-.','color',[0.5 0.5 0.5], 'linewidth',1)

axis([1 30 1e8 1e10])
xlabel('t (Days)')
ylabel('D^2 (t) (m^2)')
set(gca, 'fontsize',18)

%if j==1
%    saveas(gcf,'../figures/fit_disp_mod_shallow.eps', 'epsc')
%else
%    saveas(gcf,'../figures/fit_disp_mod_deep.eps', 'epsc')
%end


%% Plots of kurtosis with fits.

% go to plot_rel_pdf_all.m for kurtosis plots.


%% PDF plots at time when the fitting was done.
% Models
j=2; 
seps = Mod_disp(j,1).disp.^0.5;
id = fit_parms_mod(j).tstar;

maxsep = max(seps(id,:));

xax_pdf = linspace(0, maxsep,30);

[data_pdf, b]= hist(seps(id,:), xax_pdf);

norm_pdf = data_pdf./trapz(b, data_pdf);

i=3;
figure('rend','painters','pos',[10 10 800 600])
clf, hold all

bar(b/1e3, norm_pdf,'edgecolor',colors(4-j,:),'facecolor',colors(4-j,:));

ylung = lungdrendistfull(b, fit_parms_mod(j).Ro, T(id), fit_parms_mod(j).Tlung, 1);
yrich = richdistfull(b, fit_parms_mod(j).Ro, T(id), fit_parms_mod(j).beta, 1);
yray = rayleighdistfull(b, fit_parms_mod(j).Ro, T(id), fit_parms_mod(j).kappa, 1);

xlabel('$r$(km)','Interpreter','Latex')
ylabel('$2\pi r.p(r)$', 'Interpreter','latex')
box on

%axis([0 maxsep/1000 0 max([norm_pdf, yrich, ylung, yray])])
axis([0 160 0 5.5e-5])
set(gca,'FontSize',24)

legend('PDF', 'Lungdren', 'Richardson', 'Rayleigh')
legend boxoff

%if j==1
%    saveas(gcf,'../figures/pdf_model_shallow.eps', 'epsc')
%else
%    saveas(gcf,'../figures/pdf_model_deep.eps', 'epsc')
%end

%% Obs PDF
j=1

id = fit_parms_obs(j).tstar;
seps = Obs_disp(j,1).disp.^0.5;

maxsep = max(seps(id,:));

xax_pdf = linspace(0, maxsep,15);

[data_pdf, b]= hist(seps(id,:), xax_pdf);

norm_pdf = data_pdf./trapz(b, data_pdf);

i=2;
figure('rend','painters','pos',[10 10 800 600])
clf, hold all
bar(b/1e3, norm_pdf,'edgecolor',colors(4-j,:),'facecolor',colors(4-j,:));

ylung = lungdrendistfull(b, fit_parms_obs(j).Ro, T(id), fit_parms_obs(j).Tlung, 1);
yrich = richdistfull(b, fit_parms_obs(j).Ro, T(id), fit_parms_obs(j).beta, 1);
yray = rayleighdistfull(b, fit_parms_obs(j).Ro, T(id), fit_parms_obs(j).kappa, 1);

xlabel('$r$(km)','Interpreter','Latex')
ylabel('$2\pi r.p(r)$', 'Interpreter','latex')
box on
%axis([0 maxsep/1000 0 max([norm_pdf, yrich, ylung, yray])])
axis([0 160 0 5.5e-5])

set(gca,'FontSize',24)

legend('PDF', 'Lungdren', 'Richardson', 'Rayleigh')
legend boxoff

if j==1
    saveas(gcf,'../figures/pdf_obs_shallow.eps', 'epsc')
else
    saveas(gcf,'../figures/pdf_obs_deep.eps', 'epsc')
end

%% %%%%%
% PDF fits at other time


%% PDF plots at time when the fitting was done.
% Models
j=2; 
seps = Mod_disp(j,1).disp.^0.5;
id = 10;

maxsep = max(seps(id,:));

xax_pdf = linspace(0, maxsep,30);

[data_hist, bhist]= hist(seps(id,:), xax_pdf);
norm_hist = data_hist./trapz(bhist, data_hist);

[data_pdf, b] = ksdensity(seps(id,:), xax_pdf);
norm_pdf = data_pdf./trapz(b, data_pdf);
i=3;
figure('rend','painters','pos',[10 10 800 600])
clf, hold all

bar(bhist/1e3, norm_hist,'edgecolor',colors(4-j,:),'facecolor',colors(4-j,:));
%plot(b/1e3, norm_pdf,'color',colors(4-j,:), 'linewidth',4)

ylung = lungdrendistfull(b, fit_parms_mod(j).Ro, T(id), fit_parms_mod(j).Tlung, 1);
yrich = richdistfull(b, fit_parms_mod(j).Ro, T(id), fit_parms_mod(j).beta, 1);
%yray = rayleighdistfull(b, fit_parms_mod(j).Ro, T(id), fit_parms_mod(j).kappa, 1);
%yray = rayleighdistfull(b, fit_parms_mod(j).Ro, T(id), 900*24*3600, 1);

xlabel('$r$(km)','Interpreter','Latex')
ylabel('$2\pi r.p(r)$', 'Interpreter','latex')
box on

axis([0 maxsep/1500 0 max([norm_hist, yrich, ylung])])
%axis([0 160 0 5.5e-5])
set(gca,'FontSize',24)

legend('PDF', 'Lungdren', 'Richardson', 'Rayleigh')
legend boxoff

if j==1
    saveas(gcf,['../figures/pdf_model_shallow_' num2str(id) '.eps'], 'epsc')
else
    saveas(gcf,['../figures/pdf_model_deep_' num2str(id) '.eps'], 'epsc')
end

%% Obs PDF
j=2

id = 75;
seps = Obs_disp(j,1).disp.^0.5;

maxsep = max(seps(id,:));

xax_pdf = linspace(0, maxsep,30);

[data_hist, bhist]= hist(seps(id,:), xax_pdf);
norm_hist = data_hist./trapz(bhist, data_hist);


[data_pdf, b] = ksdensity(seps(id,:), xax_pdf);
norm_pdf = data_pdf./trapz(b, data_pdf);

%[data_pdf, b]= hist(seps(id,:), xax_pdf);

%norm_pdf = data_pdf./trapz(b, data_pdf);

i=2;
figure('rend','painters','pos',[10 10 800 600])
clf, hold all
bar(bhist/1e3, norm_hist,'edgecolor',colors(4-j,:),'facecolor',colors(4-j,:));
%plot(b/1e3, norm_pdf,'color',colors(4-j,:), 'linewidth',4)


ylung = lungdrendistfull(b, fit_parms_obs(j).Ro, T(id), fit_parms_obs(j).Tlung, 1);
yrich = richdistfull(b, fit_parms_obs(j).Ro, T(id), fit_parms_obs(j).beta, 1);
%yray = rayleighdistfull(b, fit_parms_obs(j).Ro, T(id), fit_parms_obs(j).kappa, 1);

xlabel('$r$(km)','Interpreter','Latex')
ylabel('$2\pi r.p(r)$', 'Interpreter','latex')
box on
axis([0 maxsep/1500 0 max([norm_hist, yrich, ylung])])
%axis([0 160 0 5.5e-5])

set(gca,'FontSize',24)

legend('PDF', 'Lungdren', 'Richardson', 'Rayleigh')
legend boxoff

if j==1
    saveas(gcf,['../figures/pdf_obs_shallow' num2str(id) '.eps'], 'epsc')
else
    saveas(gcf,['../figures/pdf_obs_deep' num2str(id) '.eps'], 'epsc')
end