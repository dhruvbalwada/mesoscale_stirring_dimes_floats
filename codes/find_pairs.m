% Last modified 29 August 2018
% Dhruv Balwada
%
% find_pairs.m
% 
% Goal of this code is to identify pairs based on some criterion of 
% separation in horizontal and vertical directions, along with ensuring
% that sufficient data is available. 

function [dispersion, correlation, pairs, pos] =find_pairs(sep, distance_class, plevel,pdiff)

ndays = 100;
comp_pres = ndays; % how many days do we compare the pressure for
separation_time = ndays;

% Save first day data (original pairs)
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
    
    % ids for pairs that are in the distance class on day 1 itself.
    id = find(origdist>=dist1 & origdist<=dist2 ...
        & meanpresdiff<=pdiff ...
        & meanorigpres>=plevel(1) ...
        & meanorigpres<=plevel(2) ...
        & (nsamples)>0.5*ndays);
    
    % ids for all other pairs, which was further away on day 1.
    idnot = find(~(origdist>=dist1 & origdist<=dist2 ...
        & meanpresdiff<=pdiff ...
        & meanorigpres>=plevel(1) ...
        & meanorigpres<=plevel(2) ...
        & (nsamples)>0.5*ndays)); % find ids of pairs that are not considered original pairs
    
    % find chance pairs
    % that are not part of original pairs
    nchanc=1;
    clear chnc_pair startid_chnc endid_chnc
    for n=1:length(idnot)
        
        % ids that are withing the pressure ranges.
        idps = find(sep(idnot(n)).dP<=pdiff  ...
                    & 0.5*(sep(idnot(n)).P1+sep(idnot(n)).P2)<=plevel(2) ...
                    &  0.5*(sep(idnot(n)).P1+sep(idnot(n)).P2)>=plevel(1));
        
        % find the minimum separation between float pairs.
        mindist = nanmin(sep(idnot(n)).dist(idps));
        
        % if minimum distance is larger than the range then move to next
        % pair.
        if ~(mindist<=dist2 && mindist>=dist1)
            continue
        end
        
        % start id corresponding to minimum distance.
        ids = find(sep(idnot(n)).dist(idps) == mindist,1);
        id1 = idps(ids);
        
        idend = find(~isnan(sep(idnot(n)).dist),1,'last');
        
        % check if sufficient samples are present in the time series of pair separations.  
        if ~isempty(id1)
            if idend-id1>ndays
                nsamples = length(find(~isnan(sep(idnot(n)).dist(id1:id1+ndays))));
            else
                nsamples = 0;
            end
            
            if nsamples>=0.5*ndays
                pres_diff = sep(idnot(n)).dP(id1);
                mean_pres_diff = nanmean(sep(idnot(n)).dP(id1:id1+comp_pres));
                pres_chan = sep(idnot(n)).P1(id1);
                if (pres_chan>=plevel(1) & pres_chan<=plevel(2) & mean_pres_diff<=pdiff)
                    chnc_pair(nchanc) = idnot(n);
                    startid_chnc(nchanc) = id1;
                    endid_chnc(nchanc) = idend;
                    nchanc = nchanc+1;
                end
            end
            
        end
    end
    
    % find pairs that have already been selected as chance or original but
    % come close again
    
    % for the floats that are original pairs
    %     for i =1:length(id)
    %         idclose = find(sep(id(i)).dist<=dist2 & sep(id(i)).dist>=dist1);
    %         idend = find(~isnan(sep(id(i)).dist),1,'last');
    %
    %         idchnc= find(idclose>=startid(id(i))+ndays,1);
    %
    %         id1 = idclose(idchnc);
    %         if idend-id1>ndays
    %             nsamples = length(find(~isnan(sep(id(i)).dist(id1:id1+ndays))));
    %         else
    %             nsamples = 0;
    %         end
    %         if nsamples>=0.75*ndays
    %             pres_diff = sep(id(i)).dP(id1);
    %             mean_pres_diff = nanmean(sep(id(i)).dP(id1:id1+comp_pres));
    %             pres_chan = sep(id(i)).P1(id1);
    %             if (pres_chan>=plevel(1) & pres_chan<=plevel(2) & mean_pres_diff<=pdiff)
    %                 chnc_pair(nchanc) = id(i);
    %                 startid_chnc(nchanc) = id1;
    %                 endid_chnc(nchanc) = idend;
    %                 nchanc = nchanc+1;
    %             end
    %         end
    %     end
    
    %     for i =1:length(idnot)
    %         idclose = find(sep(idnot(i)).dist<=dist2 & sep(idnot(i)).dist>=dist1);
    %         idend = find(~isnan(sep(idnot(i)).dist),1,'last');
    %
    %         idchnc= find(idclose>=startid(id(i))+ndays,1);
    %
    %         id1 = idclose(idchnc);
    %         if idend-id1>ndays
    %             nsamples = length(find(~isnan(sep(id(n)).dist(id1:id1+ndays))));
    %         else
    %             nsamples = 0;
    %         end
    %         if nsamples>=0.75*ndays
    %             pres_diff = sep(id(i)).dP(id1);
    %             mean_pres_diff = nanmean(sep(id(i)).dP(id1:id1+comp_pres));
    %             pres_chan = sep(id(i)).P1(id1);
    %             if (pres_chan>=plevel(1) & pres_chan<=plevel(2) & mean_pres_diff<=pdiff)
    %                 chnc_pair(nchanc) = id(n);
    %                 startid_chnc(nchanc) = id1;
    %                 endid_chnc(nchanc) = idend;
    %                 nchanc = nchanc+1;
    %             end
    %         end
    %     end
    
    % Save original pairs
    for i =1:length(id)
        sep_final(i).X(1:endid(id(i))-startid(id(i))+1) = sep(pairno(id(i))).X(startid(id(i)):endid(id(i)));
        sep_final(i).Y(1:endid(id(i))-startid(id(i))+1) = sep(pairno(id(i))).Y(startid(id(i)):endid(id(i)));
        sep_final(i).dist(1:endid(id(i))-startid(id(i))+1) = sep(pairno(id(i))).dist(startid(id(i)):endid(id(i)));
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
    
    % Save chance pairs
    n=length(id);
    for i=1:length(chnc_pair)
        sep_final(n+i).X(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).X(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).Y(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).Y(startid_chnc(i):endid_chnc(i));
        sep_final(n+i).dist(1:endid_chnc(i)-startid_chnc(i)+1) = sep(chnc_pair(i)).dist(startid_chnc(i):endid_chnc(i));
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
    
    flag = 1;
    if flag == 1
        [dispersion(l)] = rel_disp(sep_final,ndays,23);
        [correlation(l)] = correlation_func(sep_final,ndays);
        
    end
    pairs(l)= length(sep_final);
    
    
    pos(l).A = sep_final;
    %     pos(l).X2 = sep_final(:).X2;
    %     pos(l).Y1 = sep_final(:).Y1;
    %     pos(l).Y2 = sep_final(:).Y2;
    
    %     if l ==1
    %         fitter_plotter(sep_final, ndays, 23)
    %     end
    clear sep_final
    % Choose one distance class to plot the anisotropy plot for
    %     if l ==2
    %         anisotropyvst_plotter(sep_final,ndays)
    %     end
end

% relative_diff_plotter(dispersion)
