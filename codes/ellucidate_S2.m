% The goal of this appendix is to show how finite inertial ranges of energy
% spectra influence the structure functions. 

close all 
clear all 

%%
% create an idealized energy spec (piecewise 3 component)
% 
k = logspace(-2, 1, 100); 
E = zeros(size(k)); 

ka = 2*pi/200;
kb = 2*pi/15; 

% large scale constant
id = find(k<ka); 
E(id) = 1; 

% inverse energy range 

id = find(k>=ka & k<=kb); 
C1 = E(id(1)-1)/k(id(1))^(-5/3); 

E(id) = C1*k(id).^(-5/3); 

% forward enstrophy range 

id = find(k>kb); 
C2 = (C1*k(id(1)).^(-5/3))/k(id(1))^(-3);
E(id) = C2*k(id).^(-3); 

E2=E; 
C22 = (C1*k(id(1)).^(-5/3))/k(id(1))^(-4);
E2(id) = C22*k(id).^(-4); 

%%
figure('rend','painters','pos',[10 10 800 600])
g(1) = loglog(k, E, 'linewidth',3); 
hold all 
g(2) = loglog(k, E2, '--','linewidth',3); 
g(3) = loglog(k, 2*C1*k.^(-5/3),'--','color',[0.5 0.5 0.5]);
g(4) = loglog(k, 2*C2*k.^(-3),'-','color',[0.5 0.5 0.5]);

axis([1e-2 10 9e-7 1.5])
grid on 

legend(g,{'E1(k)', 'E2(k)', 'k^{-5/3}', 'k^{-3}'})
legend boxoff

set(gca,'Fontsize',22)
xlabel('k')
ylabel('E (k)')
saveas(gcf,'../figures/2D_spec.eps', 'epsc')


%% Do the calculation of the structure functions. 
r = 2*pi./k; 

for i=1:length(r) 
    S2(i) = trapz(k, E.*(1- besselj(0, r(i)*k))); 
    S22(i) = trapz(k, E2.*(1- besselj(0, r(i)*k))); 
end


%% Plot structure functions 
%close all 
figure('rend','painters','pos',[10 10 800 600])

h(1) = loglog(r, S2, 'LineWidth', 3);
grid on
hold all
h(2) = loglog(r, S22, '--', 'LineWidth', 3);

h(3) = loglog(r, 10^-(3)*r.^2, 'color', [0.5 0.5 0.5]);
%loglog(r, 10^-(2.8)*r.^1, '-.','color', [0.5 0.5 0.5])
h(4) = loglog(r, 10^-(2.35)*r.^(2/3), '--', 'color', [0.5 0.5 0.5]);

axis([min(r) max(r) min(S2) max(S2)+0.1*max(S2)])

ra = 2*pi/ka; 
loglog([ra ra] , [min(S2) max(S2)], '--', 'linewidth',1, 'color', 'r')

rb = 2*pi/kb; 
loglog([rb rb] , [min(S2) max(S2)], '--', 'linewidth',1, 'color', 'r')

set(gca, 'FontSize', 20)
xlabel('r') 
ylabel('S2 (r)')
A = legend(h, {'S1_2', 'S2_2', 'r^2', 'r^{2/3}'});
set(A, 'location', 'best')
legend boxoff
saveas(gcf,'../figures/2D_S2.eps', 'epsc')

%%
% create an idealized energy spec (piecewise 4 component)
% 
k = logspace(-2, 4, 100); 
E = zeros(size(k)); 

ka = 2*pi/200;
kb = 2*pi/15;
kc = 2*pi/1.5; 

% large scale constant
id = find(k<ka); 
E(id) = 1; 

% inverse energy range 
id = find(k>=ka & k<=kb); 
C1 = E(id(1)-1)/k(id(1))^(-5/3); 

E(id) = C1*k(id).^(-5/3); 

% forward enstrophy range 
id = find(k>kb & k<=kc); 
C2 = (C1*k(id(1)).^(-5/3))/k(id(1))^(-3);
E(id) = C2*k(id).^(-3);

% wave range 
id = find(k>kc);
C3 = (C2*k(id(1)).^(-3))/k(id(1))^(-5/3); 
E(id) = C3*k(id).^(-5/3); 

%%
%%
figure('rend','painters','pos',[10 10 800 600])
g(1) = loglog(k, E, 'linewidth',3); 
hold all 
g(2) = loglog(k, 2*C1*k.^(-5/3),'--','color',[0.5 0.5 0.5]);
g(3) = loglog(k, 2*C2*k.^(-3),'-','color',[0.5 0.5 0.5]);

axis([1e-2 100 1e-7 1.5])
grid on 

legend(g,{'E3(k)', 'k^{-5/3}', 'k^{-3}'})
legend boxoff

set(gca,'Fontsize',22)
xlabel('k')
ylabel('E (k)')
saveas(gcf,'../figures/2D_wave_spec.eps', 'epsc')

%% Do the calculation of the structure functions. 
r = 2*pi./k; 

for i=1:length(r) 
    S2(i) = trapz(k, E.*(1- besselj(0, r(i)*k))); 
end

%% Plot structure functions 
%close all 
clear h
figure('rend','painters','pos',[10 10 800 600])

h(1) = loglog(r, S2, 'LineWidth', 3);
grid on
hold all
h(2) = loglog(r, 10^-(3.2)*r.^2, '-', 'color', [0.5 0.5 0.5]);
h(3) = loglog(r, 10^-(2.3)*r.^(2/3), '--', 'color', [0.5 0.5 0.5]);

axis([min(r) max(r) min(S2) max(S2)+0.1*max(S2)])

ra = 2*pi/ka; 
loglog([ra ra] , [min(S2) max(S2)], '--', 'linewidth',1, 'color', 'r')

rb = 2*pi/kb; 
loglog([rb rb] , [min(S2) max(S2)], '--', 'linewidth',1, 'color', 'r')

rc = 2*pi/kc; 
loglog([rc rc] , [min(S2) max(S2)], '--', 'linewidth',1, 'color', 'r')

A = legend(h, {'S3_2', 'r^2', 'r^{2/3}'});
set(A, 'location', 'northwest')
legend boxoff

xlabel('r') 
ylabel('S2 (r)')

set(gca,'FontSize', 20)
saveas(gcf,'../figures/2Dwave_S2.eps', 'epsc')

close all