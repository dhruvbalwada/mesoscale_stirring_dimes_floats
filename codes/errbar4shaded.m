% first axis is time, second is bootstrap,
function [Amean, Aci] = errbar4shaded(A)
Amean = squeeze(nanmean(A,2));

for i = 1:size(A,4)
    for j =1:size(A,3)
        for k =1:size(A,1)
            Acitemp(k,:,j,i) = prctile(A(k,:,j,i), [95, 5]);
        end
    end
end
    
    for i =1:size(A,4)
        for j =1:size(A,3)
            Aci(:,1,j,i) = Acitemp(:,1,j,i) - Amean(:,j,i);
            Aci(:,2,j,i) = -Acitemp(:,2,j,i) + Amean(:,j,i);
        end
    end
end