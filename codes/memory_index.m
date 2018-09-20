% memory_index.m
%
% Function calculate the memory index
%
% Equation 5.1 in Foussard et al 2017
% M(t) = <Y.Y_0>/ Y_0<Y^2>^0.5
% Where Y is the separation vector, and Y_o is the initial separation
% vector

function [M, Mvel] = memory_index(sep,ndays)

npairs = length(sep);
X      = nan*ones(ndays,npairs);
Y      = nan*ones(ndays,npairs);
dist   = nan*ones(ndays,npairs);

for i = 1:npairs
    len = length(sep(i).X);
    if ndays<=len
        
        X(1:ndays,i) = sep(i).X(1:ndays);
        Y(1:ndays,i) = sep(i).Y(1:ndays);
        dist(1:ndays,i) = sep(i).dist(1:ndays);

        XXo(:,i) = X(:,i)*X(1,i);
        YYo(:,i) = Y(:,i)*Y(1,i);
        
        dU(1:ndays,i) = sep(i).dU(1:ndays);
        dV(1:ndays,i) = sep(i).dV(1:ndays);
        
        UUo(:,i) = dU(:,i)*dU(1,i);
        VVo(:,i) = dV(:,i)*dV(1,i);
        
    else
        
    end
end

dist0 = nanmean(sqrt(X(1,:).^2 + Y(1,:).^2));
distmean = nanmean(dist,2); 

M = nanmean(XXo + YYo, 2)/dist0./distmean; 

du0 = nanmean(dU(1,:).^2 + dV(1,:).^2).^0.5; 
dumean = nanmean(dU.^2 + dV.^2, 2).^0.5; 

Mvel = nanmean(UUo + VVo, 2)./du0./dumean; 

end
