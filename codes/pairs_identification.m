% Code to find the
%    [d,n1,nj,nk,da,n1a,n1b]=reldispbefore(Px,Py,L,latconv)
%
%    for rel. disp. calculations (to be followed by reldispafter)
%    Px,Py are the x,y matrices (each row at a given time)
%    d are the minimum distances between all pairs
%    da are the first (deployed) distances between pairs
%    the n's are indices
%    n1b is the first time when the distance falls below L(km)
%    latconv=y if want to convert to distances
%    days- the number of days for which you are calculating the dispersion


function [da,nj,nk, cnt0]=pairs_identification(Px,Py,L1,L2,latconv)

[N,M]=size(Px); % N=No. of days M=No. of floats

if(latconv=='y')
    cor=(1.1120e2)^2;                 % for converting lat/long to dist(km)
else
    cor=1;
end

% d=NaN*ones(M*(M-1)/2,1); % Min distance b/w pairs gives the option of having
% nj=d;                    % Index of first float for the minimum n1 above
% nk=d;                    % Index of second float for minimum n1 above
% da=d;                    % First deployed distance b/w pairs
cnt0=1;                  %

for j1=1:M-1            % for loops to run over all the pairs
    for j2=j1+1:M
        x1=Px(1,j1);        % X locations of floats 1
        x2=Px(1,j2);        % X locations of float with which trying to find pairings
        y1=Py(1,j1);
        y2=Py(1,j2);
        
        if(latconv=='y')    % Find the distance b/w particles on each time step
            dist=sqrt(cor)*sqrt(((x1-x2).*cosd((y1+y2)/2)).^2 + (y1-y2).^2);
        else
            dist=sqrt((x1-x2).^2 + (y1-y2).^2);
        end
                       
        if(~isempty(dist) & dist<L2 & dist>=L1) % Checks if g and dt are not empty
            
            nj(cnt0)=j1;                     % index of the first float
            nk(cnt0)=j2;                     % index of the second float
            da(cnt0)=dist;                   % dist b/w them when deployed
            cnt0=cnt0+1;                     % increase count to find next pair( keeping the same j1 go to next j2)
        end
    end
end

