% Same as load pairs, with the exception that it does not bring the first
% non-nan value to day 1. So all the particles have same time axis.

function [xp, yp,up, vp,tm, depth] = loadpairs(dep)
load ../../data/ALLDEPTHSjoinedfloats_day0to735valid
data=data0valid;

xp=[];
yp=xp;
up=xp;
vp=xp;

zf=squeeze(data(:,dep,:,1,5));
nn=min(find(isnan(zf(:,1))==0));
depth=mean(zf(nn,:))

for j=1:12
    ids = (j-1)*10+1;
    xf=squeeze(data(:,dep,:,j,3));
    yf=squeeze(data(:,dep,:,j,4));
    uf=squeeze(data(:,dep,:,j,7));
    vf=squeeze(data(:,dep,:,j,8));
    
%     xf(1:end-ids+1,:) = xf(ids:end,:);
%     xf(end-ids+2:end,:) = NaN;
%     yf(1:end-ids+1,:) = yf(ids:end,:);
%     yf(end-ids+2:end,:) = NaN;
%     uf(1:end-ids+1,:) = uf(ids:end,:);
%     uf(end-ids+2:end,:) = NaN;    
%     vf(1:end-ids+1,:) = vf(ids:end,:);
%     vf(end-ids+2:end,:) = NaN;
    
    xp=[xp xf];
    yp=[yp yf];    
    up=[up uf];
    vp=[vp vf];
end

[N,M]=size(xp);
dt=1;
tm=[0:N-1]'*dt;
end


