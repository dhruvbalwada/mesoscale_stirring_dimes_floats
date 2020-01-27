% Last Modified 29 August 2018
% Dhruv Balwada
%
% Calculate_separation_timeseries.m
% 
% This function calculates the time series of separation of all possible
% float pairs, regardless of initial separation.
% 
% It returns a structure with index number corresponding to each possible
% pair, and variables that correspond to the 2 floats that compose the
% pair.
%

function [sep] = calculate_seperation_timeseries(traj)

nflts = size(traj.Xc,2);
ndays = size(traj.Xc,1);

n=1;
% npairs = ceil(factorial(nflts)/factorial(nflts-2)/factorial(2));
npairs = nflts*(nflts-1);
sep(npairs).X = nan*ones(ndays,1);
sep(npairs).Y = nan*ones(ndays,1);
sep(npairs).dist = nan*ones(ndays,1);
sep(npairs).P1 = nan*ones(ndays,1);
sep(npairs).P2 = nan*ones(ndays,1);
sep(npairs).dP = nan*ones(ndays,1);
sep(npairs).T1 = nan*ones(ndays,1);
sep(npairs).T2 = nan*ones(ndays,1);
sep(npairs).X1 = nan*ones(ndays,1);
sep(npairs).Y1 = nan*ones(ndays,1);
sep(npairs).X2 = nan*ones(ndays,1);
sep(npairs).Y2 = nan*ones(ndays,1);
sep(npairs).U1 = nan*ones(ndays,1);
sep(npairs).V1 = nan*ones(ndays,1);
sep(npairs).U2 = nan*ones(ndays,1);
sep(npairs).V2 = nan*ones(ndays,1);
sep(npairs).dU = nan*ones(ndays,1);
sep(npairs).dV = nan*ones(ndays,1);
sep(npairs).time = nan*ones(ndays,1);
sep(npairs).names = [];

for i=1:nflts-1
    for j = i+1:nflts
        
        sep(n).X = (traj.Xc(:,i) - traj.Xc(:,j)).*cosd(0.5*(traj.Yc(:,i)+traj.Yc(:,j)))*111321;
        sep(n).Y = (traj.Yc(:,i) - traj.Yc(:,j))*111321;
        sep(n).dist = sqrt(sep(n).X.^2 + sep(n).Y.^2);
        sep(n).dU = traj.Uc(:,i) - traj.Uc(:,j); 
        sep(n).dV = traj.Vc(:,i) - traj.Vc(:,j); 
        
        sep(n).X1 = traj.Xc(:,i);
        sep(n).X2 = traj.Xc(:,j);
        sep(n).Y1 = traj.Yc(:,i);
        sep(n).Y2 = traj.Yc(:,j);
        sep(n).U1 = traj.Uc(:,i);
        sep(n).U2 = traj.Uc(:,j);
        sep(n).V1 = traj.Vc(:,i);
        sep(n).V2 = traj.Vc(:,j);
        sep(n).P1 = traj.Pi(:,i);
        sep(n).P2 = traj.Pi(:,j);
        sep(n).dP = abs(sep(n).P1 - sep(n).P2);
        sep(n).T1 = traj.Ti(:,i);
        sep(n).T2 = traj.Ti(:,j);
        sep(n).names(1) = str2num(traj.name(:,i)');
        sep(n).names(2) = str2num(traj.name(:,j)');
        sep(n).time = traj.Tv;
        n = n+1;
    end
end
