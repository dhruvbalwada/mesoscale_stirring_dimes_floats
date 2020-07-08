% Why do we need to run this convoluted code?
% It is because the output of particles is saved in a weird fashion.
% There are 12 sets of 100 particles each. Each set is released after 10
% days. So set one is released at day 1, set two at day 11 and so on.
% However, the data is saved such that all the initial positions for the
% different sets are indexed with 0. So, we need to make sure that
% particles that are present at different times are not being considered as
% members of a pair.

function [mod_sep] = idealmodel_sep_cals(mod_traj, dists)

dist1 = dists(1);
dist2 = dists(2);


tj =[];
tk =[];
tdist = [];

%for i=1:12
%    ids = (i-1)*100+1;
%    ide = i*100;
%    [dist, nj, nk, num] = pairs_identification(mod_traj.X(:,ids:ide), mod_traj.Y(:,ids:ide), dist1,dist2, 'y');
[dist, nj, nk, num] = pairs_identification(mod_traj.X, mod_traj.Y, dist1,dist2, 'n');
    %tj = [tj nj+(i-1)*100];
    
    %tk = [tk nk+(i-1)*100];
    %tdist = [tdist dist];
tj = nj;
tk = nk;
tdist = dist;
%end

%
cor=1;

for i=1:length(tdist)
    x1 = mod_traj.X(:,tj(i));
    x2 = mod_traj.X(:,tk(i));
    y1 = mod_traj.Y(:,tj(i));
    y2 = mod_traj.Y(:,tk(i));
    
    separation=sqrt((x1-x2).^2 + (y1-y2).^2);
    
    u1 = mod_traj.U(:,tj(i));
    u2 = mod_traj.U(:,tk(i));
    v1 = mod_traj.V(:,tj(i));
    v2 = mod_traj.V(:,tk(i));
    

    sep_final(i).X = (x1 - x2) ;
    sep_final(i).Y = (y1 - y2);
    sep_final(i).dist= separation;
    
    sep_final(i).dU = u1 - u2;
    sep_final(i).dV = v1 - v2;
    
    sep_final(i).X1 = x1;
    sep_final(i).Y1 = y1;
    sep_final(i).X2 = x2;
    sep_final(i).Y2 = y2;
    
    sep_final(i).U1 = u1;
    sep_final(i).V1 = v1;
    sep_final(i).U2 = u2;
    sep_final(i).V2 = v2;
    
    sep_final(i).T  = mod_traj.T;
end

mod_sep.sep = sep_final;
%mod_sep.depth = depth(dlvl);

