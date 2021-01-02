% We want to check if the form of richardson dispersion given in Graff et
% al 2015 and that proposed by Ollitrault et al 2005 look similar. 
%
% The difference is in the fact that the Ollitrault paper they assume that
% <r^2> == r^2. 
% 
% *** We find that the heurestic and more advanced curves have the same
% shape but unknown constants in front of them ...

rich_disp = @(bguess, Ro, xdata) factorial(5)/2* ...
                        ((4*bguess*xdata/9).^3) ...
                        .*exp(-9*Ro^(2/3)/4/bguess./xdata) ...
                        .*hypergeom(6,3,9*Ro^(2/3)/4/bguess./xdata);

rich_disp_naive = @(bguess, Ro, xdata) (Ro^(2/3) + 2/3*bguess*xdata).^3; 


%%
xdata = [0:100]; 
bguess =30; 
Ro =20*1e3; 

%% 
disp1 = rich_disp(bguess, Ro, xdata)/1e6; 

%%
disp2 = rich_disp_naive(bguess*3, Ro, xdata)/1e6; 

%% 

figure 
plot(xdata, disp1)
hold all 
plot(xdata, disp2)

%% 

figure 

loglog(xdata, disp1.^(1/3) - (Ro/1e3)^(2/3), 'linewidth',2)
hold all 
loglog(xdata, disp2.^(1/3) - (Ro/1e3)^(2/3), 'linewidth',2)
loglog(xdata, xdata,'--')
legend('Exact', 'Ollitrault', 'x=x')
xlabel('Days')
set(gca,'fontsize',18)

ylabel('[$\overline{r^{2}}^{1/3}(t | r_0) - r_0^{2/3}] (km^{2/3})$', 'Interpreter','Latex')
title('Ro = 20km, \beta=30m^{2/3}/day')
