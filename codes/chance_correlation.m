% chance_correlation.m
% 
% Function calculates the correlation of initial separation to initial
% separation velocity.
% 
% The term to be calculated is <dV.Y_o>
%
% Hypothesis 1 of Babiano et al 1990 assumes that this is 0, which is the
% case of the pairs are all initial randomly oriented regardless of 
% velocity field (original pairs).

function [C] = chance_correlation(sep, ndays)

npairs = length(sep);
X      = nan*ones(ndays,npairs);
Y      = nan*ones(ndays,npairs);
dist   = nan*ones(ndays,npairs);

for i = 1:npairs
    len = length(sep(i).X);
    if ndays<=len
        
        X(1:ndays,i) = sep(i).X(1:ndays);
        Y(1:ndays,i) = sep(i).Y(1:ndays);
        
        dU(1:ndays,i) = sep(i).dU(1:ndays);
        dV(1:ndays,i) = sep(i).dV(1:ndays);

        dUXo(:,i) = dU(:,i)*X(1,i);
        dVYo(:,i) = dV(:,i)*Y(1,i);
        
    else
        
    end
end

dist0 = nanmean(sqrt(X(1,:).^2 + Y(1,:).^2));
dvelmean = nanmean(dU.^2 + dV.^2, 2);

C = nanmean(dUXo + dVYo, 2)/dist0./sqrt(dvelmean);

end