function [ specmean, specci ] = spectra_logerrors( spectra )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
nsamps = size(spectra, 2); 
nfreq = size(spectra, 1); 

for n =1:500
    y = randsample(nsamps, nsamps, true); 
   
    specbs(:,n) = nanmean(spectra(:,y), 2); 
end

specmean = nanmean(specbs, 2); 

for i = 1:nfreq
    spectemp(i,:) = prctile(specbs(i,:),[95, 5]); 
end

specci(:,1) = 0.434*(spectemp(:,1) - specmean)./specmean; 
specci(:,2) = 0.434*(-spectemp(:,2) + specmean)./specmean; 

end

