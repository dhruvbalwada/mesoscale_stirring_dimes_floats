% Last modified 29 August 2018
% Dhruv Balwada
%
% find_pairs.m
%
% Goal of this code is to identify pairs based on some criterion of
% separation in horizontal and vertical directions, along with ensuring
% that sufficient data is available.

function [pairs, pos] =find_pairs(sep, distance_class, plevel,pdiff, ndays)

comp_pres = ndays; % how many days do we compare the pressure for
separation_time = ndays;

for i =1:length(sep)
    id = find(~isnan(sep(i).dist),1);
    
    if ~isempty(id)
        origdist(i)     = sep(i).dist(id);
        meanorigpres(i) = 0.5*(sep(i).P1(id)+sep(i).P2(id));
        origpresdiff(i) = sep(i).dP(id);
        meanpresdiff(i) = nanmean(sep(i).dP(id:id+comp_pres));
        pairno(i)       = i;
        startid(i)      = id;
        endid(i)        = find(~isnan(sep(i).dist),1,'last');
        nsamples(i)     = length(find(~isnan(sep(i).dist(startid(i):startid(i)+ndays))));
    else
        origdist(i)     = NaN;
        pairno(i)       = NaN;
        startid(i)      = NaN;
        endid(i)        = NaN;
        meanorigpres(i) = NaN;
        origpresdiff(i) = NaN;
        meanpresdiff(i) = NaN;
        nsamples(i)     = NaN;
    end
    
end


% Start adding pairs to different distance classes
for l=1:length(distance_class)-1
    
    dist1 = distance_class(l);     % D_o class
    dist2 = distance_class(l+1);     % D_o
    
    % ids for all pairs, which satisfy pressure criterion
    % do a basic check first.
    idnot = find((meanpresdiff <= pdiff ...
        & meanorigpres   >= plevel(1) ...
        & meanorigpres   <= plevel(2) ...
        & nsamples       >  0.25*ndays));
    
    % find pairs
    % that are not part of original pairs
    nchanc=1;
    clear chnc_pair startid_chnc endid_chnc
    
    for n=1:length(idnot)
        
        % ids that are within the pressure ranges.
        idps = find(sep(idnot(n)).dP<=pdiff  ...
            & 0.5*(sep(idnot(n)).P1+sep(idnot(n)).P2)<=plevel(2) ...
            & 0.5*(sep(idnot(n)).P1+sep(idnot(n)).P2)>=plevel(1));
        
        % find the minimum separation between float pairs.
        mindist = nanmin(sep(idnot(n)).dist(idps));
        
        % if minimum distance is larger than the range then move to next
        % pair.
        
        %if ~(mindist<=dist2 & mindist>=dist1)
        %    continue
        %end
        
        % start id corresponding to first time distance comes in range
        
        ids = find(sep(idnot(n)).dist(idps) >=dist1 & ...
            sep(idnot(n)).dist(idps) <= dist2 , 1);
        
        id1 = idps(ids);
        
        idend = find(~isnan(sep(idnot(n)).dist),1,'last');
        
        % check if sufficient samples are present in the time series of pair separations.
        if ~isempty(id1)
            if idend-id1>ndays
                nsamples = length(find(~isnan(sep(idnot(n)).dist(id1:id1+ndays))));
            else
                nsamples = 0;
            end
            
            if nsamples>=0.25*ndays
                
                pres_diff      = sep(idnot(n)).dP(id1);
                mean_pres_diff = nanmean(sep(idnot(n)).dP(id1:id1+comp_pres));
                pres_chan      = sep(idnot(n)).P1(id1);
                
                if (pres_chan>=plevel(1) & pres_chan<=plevel(2) & mean_pres_diff<=pdiff)
                    chnc_pair(nchanc)    = idnot(n);
                    startid_chnc(nchanc) = id1;
                    endid_chnc(nchanc)   = idend;
                    
                    nchanc = nchanc+1;
                end
            end
            
        end
    end
    

    flag =1;
    
    if flag ==1
        chnc_pair_temp = chnc_pair;
        
        for i =1:length(chnc_pair_temp) % cycle through all chance pairs
            
            % find ids when the distance bw pairs is within criterion
            idclose = find(sep(chnc_pair_temp(i)).dist<=dist2 & sep(chnc_pair_temp(i)).dist>=dist1);
            idend = find(~isnan(sep(chnc_pair_temp(i)).dist), 1,'last');
            
            % find an id that is away from the chance pair id.
            idchnc= find(idclose>=startid_chnc(i)+0.25*ndays, 1);
            
            if isempty(idchnc)
                continue
            end
            
            id1 = idclose(idchnc);
            
            if idend-id1>ndays
                nsamples = length(find(~isnan(sep(chnc_pair_temp(i)).dist(id1:id1+ndays))));
            else
                nsamples = 0;
            end
            
            if nsamples>=0.25*ndays
                pres_diff = sep(chnc_pair_temp(i)).dP(id1);
                mean_pres_diff = nanmean(sep(chnc_pair_temp(i)).dP(id1:id1+comp_pres));
                pres_chan = sep(chnc_pair_temp(i)).P1(id1);
                if (pres_chan>=plevel(1) & pres_chan<=plevel(2) & mean_pres_diff<=pdiff)
                    chnc_pair(nchanc) = chnc_pair_temp(i);
                    startid_chnc(nchanc) = id1;
                    endid_chnc(nchanc) = idend;
                    nchanc = nchanc+1;
                end
            end
        end
    end

    % Save chance pairs
    n=0;
    
    for i=1:length(chnc_pair)
        sep_final(n+i).X(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).X(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).Y(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).Y(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).dist(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).dist(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).dU(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).dU(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).dV(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).dV(startid_chnc(i):endid_chnc(i));
        
        sep_final(n+i).U1(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).U1(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).U2(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).U2(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).V1(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).V1(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).V2(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).V2(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).X1(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).X1(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).X2(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).X2(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).Y1(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).Y1(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).Y2(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).Y2(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).P1(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).P1(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).P2(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).P2(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).names = sep(pairno(chnc_pair(i))).names;
        
    end
        
    pairs(l)= length(sep_final);
    
    pos(l).sep = sep_final;
      
end
