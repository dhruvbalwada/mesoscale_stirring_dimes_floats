function y=rayleighdistfull(x,x0,t,K2,flag)

%dx=mean(diff(x))/10;
%x2=[0:dx:max(x)];

x2 = x;

y= besseli(0,x0*x2/2/K2/t).*exp(-(x0^2+x2.^2)/4/K2/t)./4/pi/K2/t;

% trap(y,dx);
if flag==1
plot(x2/1000,y,'-.','color',[0.5 0.5 0.5],'linewidth',2)
end
% x2 = x;
% y= 2*pi*besseli(0,x0*x2/2/K2/t).*exp(-(x0^2+x2.^2)/4/K2/t).*x2/4/pi/K2/t;
