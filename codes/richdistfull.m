
% the richardson distribution, full (non-asymptotic) version
% the separations, x, are uniformly separated
% bet is the amplitude of the diffusivity (k = bet* r^4/3)

function y=richdistfull(x,x0,t,beta, flag)

y = 3/4/pi/beta/t/x0^(2/3)./x.^(2/3).* ...
    exp(-9*(x0^(2/3)+x.^(2/3))/4/beta/t).*  ...
    besseli(2,9*x0^(1/3)*x.^(1/3)/2/beta/t);

%dx=mean(diff(x))/10;
%x2=[0:dx:max(x)];
% x2 = x;
%y=x2.^(1/3).*exp(-9*(x0^(2/3)+x2.^(2/3))/4/bet/t)...   % this is r p(r)
%    .*besseli(2,9*x0^(1/3)*x2.^(1/3)/2/bet/t);
%norm=(2*bet*t*x0^(2/3))/3;
%    norm=trap(y,dx);
%y=y/norm;

%y2=x2.*exp(-9*x2.^(2/3)/4/bet/t);
%norm=3*(4*bet*t/9)^3;
%y2=y2/norm;

% plot(x2,y,x2,y2,'--')
if flag==1
    plot(x2/1000,y,'--','color',[0.5 0.5 0.5],'linewidth',2)
end
% legend('full','asymptotic')
% legend('full')
% pause

% disp(['rich norm= ',num2str(trap(y,dx))])

% x2 = x;
% y=x2.^(1/3).*exp(-9*(x0^(2/3)+x2.^(2/3))/4/bet/t)...   % this is r p(r)
%     .*besseli(2,9*x0^(1/3)*x2.^(1/3)/2/bet/t);
% norm=(2*bet*t*x0^(2/3))/3;
% %    norm=trap(y,dx);
% y=y/norm;
