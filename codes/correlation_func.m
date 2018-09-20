function [correlation] = correlation_func(sep,ndays)

%%  Correlations in time

npairs = length(sep);
U1 = nan*ones(ndays,npairs);
V1 = nan*ones(ndays,npairs);
U2 = nan*ones(ndays,npairs);
V2 = nan*ones(ndays,npairs);
dist = nan*ones(ndays,npairs);

for i = 1:npairs
    len = length(sep(i).U1);
    if ndays<=len
        U1(1:ndays,i) = sep(i).U1(1:ndays);
        V1(1:ndays,i) = sep(i).V1(1:ndays);
        U2(1:ndays,i) = sep(i).U2(1:ndays);
        V2(1:ndays,i) = sep(i).V2(1:ndays);
        dist(1:ndays,i) = sep(i).dist(1:ndays);
    else
    end
end

for i=1:ndays
    
    cor(i) = nanmean( U1(i,:).*U2(i,:) + V1(i,:).*V2(i,:))/ ...
        nanmean(U1(i,:).^2 + V1(i,:).^2).^0.5/nanmean(U2(i,:).^2 + V2(i,:).^2).^0.5;
    
    coru(i) = nanmean(U1(i,:).*U2(i,:))./(nanmean(U1(i,:).^2).^0.5* nanmean(U2(i,:).^2).^0.5);
    corv(i) = nanmean(V1(i,:).*V2(i,:))./(nanmean(V1(i,:).^2).^0.5* nanmean(V2(i,:).^2).^0.5);
    
    r(i) = nanmean(dist(i,:));
end

correlation.c = cor;
correlation.cuu = coru;
correlation.cvv = corv;

correlation.r   = r;