% Dhruv Balwada
% 10 Jan 2020

clear all 
close all

% how do different PDFs evolve in time 

%% Diffusive 

dr = 100;
r = 0+dr:dr:300e3; 
t = (1:200)*24*3600;

beta = 0.05*(1000)^(2/3)/24/3600; % beta
r0 = 1e3; 

Pdiff = zeros(length(r), length(t)); 

for i = 1:length(t)
    Pdiff(:, i) = richdistfull(r, r0, t(i), beta, 0);
end

%% r2 second moment 
mom0 = 2*pi*sum((r.^1).*Pdiff'*dr,2);
mom1 = 2*pi*sum((r.^2).*Pdiff'*dr,2);
mom2 = 2*pi*sum((r.^3).*Pdiff'*dr,2);
mom3 = 2*pi*sum((r.^4).*Pdiff'*dr,2);
mom4 = 2*pi*sum((r.^5).*Pdiff'*dr,2);

%% 
figure
hold all 
for i =1:10:100
    plot(r, Pdiff(:,i)) 
end


%% 
figure 
subplot(121)
%hold all 
contourf(r,t, Pdiff')

%%
figure
contourf(r,t, 2*pi*cumsum(r.*Pdiff',2)*dr, 'edgecolor','none')
hold all 
plot(mom1, t, '--', 'linewidth', 2, 'color', 'k')
plot(mom2.^0.5, t, 'linewidth', 2, 'color', 'k')
xlim([0,100e3])
 
%% 
figure, 
loglog(t, mom1)
hold all 
loglog(t, mom2.^(1/2))
loglog(t, mom3.^(1/3))
loglog(t, mom4.^(1/4))

%% 
close all 
figure
plot(t, mom4./(mom2.^2))
