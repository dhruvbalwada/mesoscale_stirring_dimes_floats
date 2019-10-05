% How similar do exponential and ballistic growth curves look. 

T = [0:40] ;

Do =1; 
Tau = 6; 

Dexp = Do*exp(T/Tau); 

C=0.05; 
Dbal = Do + C*T.^2; 

%% loglog plot
figure 
loglog(T, Dexp)
hold all 
loglog(T, Dbal)

%% semilog plot
figure
semilogy(T, Dexp)
hold all 
semilogy(T, Dbal)

%% compensated 

figure
semilogx(T, (Dexp - Dexp(1))./T.^2)
hold all
semilogx(T, (Dbal - Dbal(1))./T.^2)

%% characteristic curve 
diffbal= 0*Dbal; 
diffexp= 0*Dexp; 

diffbal(2:end-1) = (Dbal(3:end) - Dbal(1:end-2))./2;
diffexp(2:end-1) = (Dexp(3:end) - Dexp(1:end-2))./2;

Taubal = Dbal./diffbal; 
Tauexp = Dexp./diffexp; 

%% 
figure 
semilogx(Dbal.^0.5, Taubal, 'o-')
hold all 
semilogx(Dexp.^0.5, Tauexp, 'o-')

%% Plot diffusivity 

figure 
loglog(Dbal.^0.5, diffbal); 
hold all 
loglog(Dexp.^0.5, diffexp); 
dist =[1,100]; 
loglog(dist, 0.2*dist.^2,'--')
