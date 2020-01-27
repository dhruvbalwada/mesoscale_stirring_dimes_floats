
% the lungdren distribution, full (non-asymptotic) version
% the separations, x, are uniformly separated
% tlung is the time scale that goes into the distribution.

function y=lungdrendistfull(x,x0,t,tlung, flag)

y = 1/4/pi/(pi*t/tlung)^(1/2)/x0^2* ... 
    exp(-(log(x./x0)+2*t/tlung).^2/4/t*tlung);

%dx=mean(diff(x))/10;
%x2=[0:dx:max(x)];
% x2 = x;
%y=x2.*exp(-(log(x2./x0)+2*t/tlung).^2/4/t*tlung);   % this is r p(r)
%norm=4*pi^1.5*(t/tlung)^0.5*x0^2/2/pi;
%    norm=trap(y,dx);
%y=y/norm;

if flag ==1
plot(x2/1000,y,'-','color',[0.5 0.5 0.5],'linewidth',2)
end
% legend('lungdren')
% pause

% x2 = x;
% y=x2.*exp(-(log(x2./x0)+2*t/tlung).^2/4/t*tlung);   % this is r p(r)
% norm=4*pi^1.5*(t/tlung)^0.5*x0^2/2/pi;
% %    norm=trap(y,dx);
% y=y/norm;

% disp(['lungdren norm= ',num2str(trap(y,dx))])
