% Last modified 29 August 2018
% Dhruv Balwada
% 
% find_pairs_random.m
%
% This function takes all possible float pairs and identifies the ones that
% are based on the criterion that are provided.


% traj - structure with the data loaded into it
% After every integral time(chosen to be 10 days here) the float pair is
% reinitialized as a new pair.


function [pairs,pos] = find_pairs_random(sep, distance_class, plevel, pdiff, ndays, int_time)

comp_pres = ndays; % how many days do we compare the pressure for

norandompairs =1;

for i =1:length(sep)
    
    id1 = find(~isnan(sep(i).dist),1); % first non-nan
    idlast = find(~isnan(sep(i).dist),1,'last'); % last non-nan.
    
    n= floor((idlast-id1 - comp_pres)/int_time); % possible number of pairs, 
                                                 % considering that a pair loses memory over ndays.
    
    if n==0, continue, end
    
    for k =1:n
        id = id1 + (k-1)*int_time;
        
        pairno(norandompairs) = i; % pair number
        startid(norandompairs) = id; % starting id of first non-nan
        endid(norandompairs) = find(~isnan(sep(i).dist),1,'last'); % ending id of last non-nan
        nsamples(norandompairs) = length(find(~isnan(sep(i).dist( ...
             startid(norandompairs):startid(norandompairs)+ndays)))); 
                                     % how many non-nans (samples) exist between start and end id.
        
        origdist(norandompairs) = sep(i).dist(id);
        meanorigpres(norandompairs) = 0.5*(sep(i).P1(id)+sep(i).P2(id));
        origpresdiff(norandompairs) = sep(i).dP(id);
        meanpresdiff(norandompairs) = nanmean(sep(i).dP(id:id+comp_pres));
        meanpres(norandompairs) = nanmean((sep(i).P1(id:id+comp_pres)+sep(i).P2(id:id+comp_pres))/2);
        norandompairs = norandompairs + 1;
        
    end
    
end


for l=1:length(distance_class)-1
    
    dist1 = distance_class(l);     % D_o class
    dist2 = distance_class(l+1);     % D_o
    
    id = find(origdist>=dist1 & origdist<=dist2 & meanpresdiff<=pdiff & meanpres>=plevel(1) ...
        & meanpres<=plevel(2) & nsamples>=0.5*ndays);
    
    pairs(l)= length(id);
    
    for i =1:length(id)
        sep_final(i).X(1:endid(id(i))-startid(id(i))+1)    = sep(pairno(id(i))).X(startid(id(i)):endid(id(i)));
        sep_final(i).Y(1:endid(id(i))-startid(id(i))+1)    = sep(pairno(id(i))).Y(startid(id(i)):endid(id(i)));
        sep_final(i).dist(1:endid(id(i))-startid(id(i))+1) = sep(pairno(id(i))).dist(startid(id(i)):endid(id(i)));
        sep_final(i).dU(1:endid(id(i))-startid(id(i))+1)   = sep(pairno(id(i))).dU(startid(id(i)):endid(id(i)));
        sep_final(i).dV(1:endid(id(i))-startid(id(i))+1)   = sep(pairno(id(i))).dV(startid(id(i)):endid(id(i)));
      
        sep_final(i).U1(1:endid(id(i))-startid(id(i))+1) = sep(pairno(id(i))).U1(startid(id(i)):endid(id(i)));
        sep_final(i).U2(1:endid(id(i))-startid(id(i))+1) = sep(pairno(id(i))).U2(startid(id(i)):endid(id(i)));
        sep_final(i).V1(1:endid(id(i))-startid(id(i))+1) = sep(pairno(id(i))).V1(startid(id(i)):endid(id(i)));
        sep_final(i).V2(1:endid(id(i))-startid(id(i))+1) = sep(pairno(id(i))).V2(startid(id(i)):endid(id(i)));
        sep_final(i).X1(1:endid(id(i))-startid(id(i))+1) = sep(pairno(id(i))).X1(startid(id(i)):endid(id(i)));
        sep_final(i).X2(1:endid(id(i))-startid(id(i))+1) = sep(pairno(id(i))).X2(startid(id(i)):endid(id(i)));
        sep_final(i).Y1(1:endid(id(i))-startid(id(i))+1) = sep(pairno(id(i))).Y1(startid(id(i)):endid(id(i)));
        sep_final(i).Y2(1:endid(id(i))-startid(id(i))+1) = sep(pairno(id(i))).Y2(startid(id(i)):endid(id(i)));        
        sep_final(i).P1(1:endid(id(i))-startid(id(i))+1) = sep(pairno(id(i))).P1(startid(id(i)):endid(id(i)));
        sep_final(i).P2(1:endid(id(i))-startid(id(i))+1) = sep(pairno(id(i))).P2(startid(id(i)):endid(id(i)));
        sep_final(i).names = sep(pairno(id(i))).names;
    end
    
     pos(l).sep = sep_final;
     
  
end