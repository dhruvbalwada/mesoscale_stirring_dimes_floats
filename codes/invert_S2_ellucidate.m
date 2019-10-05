% Here we will use an inverse method to invert for the optimal Energy. 

%% Create an artificial E 
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

close all 
figure 
subplot(121) 
plot(k, E, 'linewidth', 2)
xlabel('k') 
ylabel('E')
set(gca,'fontsize',20)

subplot(122)
loglog(k, E, 'linewidth',2)
xlabel('k') 
ylabel('E')
title('log-log')
set(gca,'fontsize',20)
%% create the corresponding S2 
% Do the calculation of the structure functions. 
r = 2*pi./k; % r is the wavelength corresponding to the wavenumber k 

S2 = S2_from_E(E, k, r); % Evaluate the an integral, S2 is an integral transform of E. 
%
figure 
subplot(121)
plot(r, S2, 'linewidth',2)
xlabel('r') 
ylabel('S2')
set(gca,'fontsize',20)
axis([min(r) max(r) min(S2) max(S2)])

subplot(122)
loglog(r, S2, 'linewidth',2)
xlabel('r') 
ylabel('S2')
title('log-log')
set(gca,'fontsize',20)
axis([min(r) max(r) min(S2) max(S2)])

%% invert S2 using the inverse relationship
% this is the inverse transform relationship

clear E_inverted

for i = 1:length(k) 
    E_inverted(i) = trapz(r, (S2(end) - S2)'/4.*besselj(0, k(i)*r).*k(i).*r); 
end
%
%% this method just gives

close all 
figure 
subplot(121) 
plot(k, E, 'linewidth', 2)
hold all 
plot(k, E_inverted, 'linewidth', 2)
xlabel('k') 
ylabel('E')
set(gca,'fontsize',20)

subplot(122)
loglog(k, E, 'linewidth',2)
xlabel('k') 
ylabel('E')
title('log-log')
set(gca,'fontsize',20)

%% If we directly pass the function that calculates S2 to PCG
E_inv = pcg(@(Ein)S2_from_E(Ein, k', r), (S2));

%%
E_inv = pcg(@(Ein)norm_S2(Ein, k', r, 1), S2, 1e-3, 1500); 

%% 
figure
subplot(121)
plot(k, E) 
hold all 
plot(k, E_inv)

subplot(122)
loglog(k, E) 
hold all 
loglog(k, E_inv)

%%
figure 
loglog(r, S2)
hold all 
loglog(r, S2_from_E(E_inv, k,r))
