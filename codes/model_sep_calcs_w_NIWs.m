% Why do we need to run this convoluted code?
% It is because the output of particles is saved in a weird fashion.
% There are 12 sets of 100 particles each. Each set is released after 10
% days. So set one is released at day 1, set two at day 11 and so on.
% However, the data is saved such that all the initial positions for the
% different sets are indexed with 0. So, we need to make sure that
% particles that are present at different times are not being considered as
% members of a pair.

function [mod_sep] = model_sep_calcs_w_NIWs(mod_traj, dists,f, A, B, C, IO_flag,L)

dist1 = dists(1)/1e3;
dist2 = dists(2)/1e3;


tj =[];
tk =[];
tdist = [];

for i=1:12
    ids = (i-1)*100+1;
    ide = i*100;
    [dist, nj, nk, num] = pairs_identification(mod_traj.X(:,ids:ide), ...
        mod_traj.Y(:,ids:ide), dist1,dist2, 'y');
    tj = [tj nj+(i-1)*100];
    tk = [tk nk+(i-1)*100];
    tdist = [tdist dist];
end

%
cor=(1.1120e5)^2;

for i=1:length(tdist)
    x1 = mod_traj.X(:,tj(i));
    x2 = mod_traj.X(:,tk(i));
    y1 = mod_traj.Y(:,tj(i));
    y2 = mod_traj.Y(:,tk(i));
    
    time = mod_traj.T;
    
    
    Arand = 2*A*rand(1);
    Brand = 2*B*rand(1);
    
    
    if IO_flag==1
        
        % Code for IOs where the IOs lose memory in time
        
        %Arand = A;
        %Brand = B;
        
        x1 = x1 + Arand*sin(f*time)./cosd((y1+y2)/2);
        x2 = x2 + ( Arand + Brand*( 1-exp(-C*time)) ).*sin(f*time)./cosd((y1+y2)/2);
        
        y1 = y1 + Arand*(cos(f*time) -1);
        y2 = y2 + ( Arand + Brand*(1-exp(-C*time)) ).*(cos(f*time) -1);
    
    elseif IO_flag ==2
    
        % Code for IOs where the IOs lose memory in space
        separation=sqrt((cor)*(((x1-x2).*cosd((y1+y2)/2)).^2 + (y1-y2).^2));
        
        rand_phase = rand(1)*2*pi;
        
        x1 = x1 + Arand*sin(f*time + rand_phase)./cosd((y1+y2)/2);
        y1 = y1 + Arand*(cos(f*time + rand_phase) -1);
        
         
        
        id1 = find(separation<L);
        x2(id1) = x2(id1) + ( Arand + Brand*((separation(id1)).^(C)) ).*sin(f*time(id1) + rand_phase)./cosd((y1(id1)+y2(id1))/2);
        y2(id1) = y2(id1) + ( Arand + Brand*((separation(id1)).^(C)) ).*(cos(f*time(id1) + rand_phase) -1);
        id2 = find(separation>=L);
        x2(id2) = x2(id2) + ( Arand + Brand*((L).^(C)) ).*sin(f*time(id2) + rand_phase)./cosd((y1(id2)+y2(id2))/2);
        y2(id2) = y2(id2) + ( Arand + Brand*((L).^(C)) ).*(cos(f*time(id2) + rand_phase) -1 );        
        
    elseif IO_flag == 3 
        % Code for IOs where the IOs lose memory in space
        separation=sqrt((cor)*(((x1-x2).*cosd((y1+y2)/2)).^2 + (y1-y2).^2));
        
        x1 = x1 + Arand*sin(f*time)./cosd((y1+y2)/2);
        y1 = y1 + Arand*(cos(f*time) -1);
        
        id1 = find(separation<L);
        x2(id1) = x2(id1) + ( Arand + Brand*((separation(id1)).^(C)) ).*sin(f*time(id1))./cosd((y1(id1)+y2(id1))/2);
        y2(id1) = y2(id1) + ( Arand + Brand*((separation(id1)).^(C)) ).*(cos(f*time(id1)) -1);
        id2 = find(separation>=L);
        %x2(id2) = x2(id2) + ( Arand + 0*Brand*((L).^(C)) ).*sin(f*time(id2))./cosd((y1(id2)+y2(id2))/2);
        %y2(id2) = y2(id2) + ( Arand + 0*Brand*((L).^(C)) ).*(cos(f*time(id2)) -1 );
        nslope = 2; 
        const = exp(log(Brand) + C*log(L) + nslope*log(L)); 
        
        x2(id2) = x2(id2) + ( Arand + const*((separation(id2)).^(-nslope) ) ).*sin(f*time(id2))./cosd((y1(id2)+y2(id2))/2);
        y2(id2) = y2(id2) + ( Arand + const*((separation(id2)).^(-nslope) ) ).*(cos(f*time(id2)) -1 );
    end
    
    separation=sqrt((cor)*(((x1-x2).*cosd((y1+y2)/2)).^2 + (y1-y2).^2));
    
    
    %u1 = mod_traj.U(:,tj(i));
    %u2 = mod_traj.U(:,tk(i));
    %v1 = mod_traj.V(:,tj(i));
    %v2 = mod_traj.V(:,tk(i));
    
    
    sep_final(i).X = (x1 - x2).*cosd((y1+y2)/2)*sqrt(cor) ;
    sep_final(i).Y = (y1 - y2)*sqrt(cor);
    sep_final(i).dist= separation;
    
    %sep_final(i).dU = u1 - u2;
    %sep_final(i).dV = v1 - v2;
    
    sep_final(i).X1 = x1;
    sep_final(i).Y1 = y1;
    sep_final(i).X2 = x2;
    sep_final(i).Y2 = y2;
    
    sep_final(i).T  = time;
    
    % get the velocity from forward differencing (will be important for structure functions)
    U1 = nan*sep_final(i).X1;
    V1 = nan*sep_final(i).X1;
    U2 = nan*sep_final(i).X1;
    V2 = nan*sep_final(i).X1;
    dt = (sep_final(i).T(2)- sep_final(i).T(1))*24*3600; 
    
    %cor = 110*1e3; 

    U1(1:end-1,:) = diff(sep_final(i).X1)*sqrt(cor).*cosd(sep_final(i).Y1(1:end-1,:))/dt; 
    V1(1:end-1,:) = diff(sep_final(i).Y1)*sqrt(cor)/dt; 
    U2(1:end-1,:) = diff(sep_final(i).X2)*sqrt(cor).*cosd(sep_final(i).Y2(1:end-1,:))/dt; 
    V2(1:end-1,:) = diff(sep_final(i).Y2)*sqrt(cor)/dt; 
    
    sep_final(i).U1 = U1;
    sep_final(i).V1 = V1;
    sep_final(i).U2 = U2;
    sep_final(i).V2 = V2;
    
    sep_final(i).dU = U1 - U2; 
    sep_final(i).dV = V1 - V2; 
end

mod_sep.sep = sep_final;
%mod_sep.depth = depth(dlvl);

