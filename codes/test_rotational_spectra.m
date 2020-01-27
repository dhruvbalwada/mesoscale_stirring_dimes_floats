clear all 
close all 
%%

A = 1; 
T = 1:1/24:100; 
f = 1e-4; 

X = A*sin(f*T*24*3600);
Y = A*(cos(f*T*24*3600)-1); 

[Xl, Yl] = xy2latlon(X,Y, -55, 200); 

%Y = -55 + A/110*(cos(f*T*24*3600)-1); 
%X = 200 + A/110*sin(f*T*24*3600); 

%% 

figure
plot(T, X)
%%
figure
plot(X,Y)

%% 

CV = latlon2uv(T, Y', X')/100; 

psi = sleptap(length(CV)); 

[F, SPP, SNN, SPN] = mspec(CV, psi);

%%
figure
plot(T, real(CV))
hold all
plot(T, imag(CV))
%% 
figure
quiver(T, 0*T, real(CV)', imag(CV)') 

%%

figure 
loglog(F, SPP)
hold all
loglog(F, SNN)