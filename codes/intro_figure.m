% intro_figure.m
% Modified by Dhruv Balwada
% on 1 September 2018

% load topography

topo=load('../../data/topo2_SS_SO.mat');
topo.lon=topo.lon-360;
topo.bath=nanmoving_average2(topo.bath,4,2);
topo.bath(topo.bath<-7000)=-7000;
loni=topo.lon(1:5:end);
lati=topo.lat(1:2:end);
[X2 Y2] = meshgrid(loni,lati);
[X1 Y1] = meshgrid(topo.lon,topo.lat);
% change the topography to smoothed version (only for speed up of ploting
% as it is just a qualitative task)
bath = interp2(X1, Y1, topo.bath', X2, Y2,'cubic');
topo_intv=[-7000:500:0];

%% 
% load Orsi Fronts
pf= load('../../data/fronts/pf.txt');
saf = load('../../data/fronts/saf.txt');

%% 
% float trajectories 

fh = figure()

set(fh, 'Visible','off');
clf
m_proj('Lambert','lat',[-66 -51], 'lon',[-114 -72]) 
days_plt = 100; 
hold all 
col = colormap(1-gray(100));

colors = lines(5); 

%m_elev('contourf', [-6000 -4000 -2000],'edgecolor','none')

m_contourf(loni, lati, bath, 'edgecolor','none')
%colorbar()

m_coast('patch',[0 0 0])
m_plot(pf(:,1), pf(:,2),'--','linewidth',1,'color',colors(4,:));
m_plot(saf(:,1), saf(:,2),'--','linewidth',1,'color',colors(4,:));

for j=1:size(traj.Xc, 2)
    ids = find(~isnan(traj.Xc(:,j)),1);
    
    mp = nanmean(traj.Pi( ids: ids+days_plt, j));
    
    
    if mp>1000
        m_plot(traj.Xc(ids: ids+days_plt, j), traj.Yc(ids:ids+days_plt, j), 'linewidth',1, 'color',colors(1,:));
    elseif mp <1000
        m_plot(traj.Xc(ids: ids+days_plt, j), traj.Yc(ids:ids+days_plt, j), 'linewidth',1, 'color',colors(2,:));
    end
    
    m_plot(traj.Xc(ids, j), traj.Yc(ids, j), '.', 'Markersize',10,'color',colors(5,:));

end 
set(gca, 'fontsize', 16)


m_grid('linest', 'none')

saveas(gca,'../figures/intro.eps')

%print(gcf, '-r0', 'intro.pdf');

%%
% Model particles

d = [4, 10];

[mod_traj1.X, mod_traj1.Y, mod_traj1.U, mod_traj1.V, mod_traj1.Z, mod_traj1.T, depth(1)] = loadpairs_with_depth(d(1));
[mod_traj2.X, mod_traj2.Y, mod_traj2.U, mod_traj2.V, mod_traj2.Z, mod_traj2.T, depth(2)] = loadpairs_with_depth(d(2));


%%
fh = figure()

set(fh, 'Visible','off');
clf
m_proj('Lambert','lat',[-66 -51], 'lon',[-114 -72]) 
days_plt = 100; 
hold all 
col = colormap(1-gray(100));

colors = lines(5); 

%m_elev('contourf', [-6000 -4000 -2000],'edgecolor','none')

m_contourf(loni, lati, bath, 'edgecolor','none')
%colorbar()

m_coast('patch',[0 0 0])
m_plot(pf(:,1), pf(:,2),'--','linewidth',1,'color',colors(4,:));
m_plot(saf(:,1), saf(:,2),'--','linewidth',1,'color',colors(4,:));

for j=1:4:240
    ids = find(~isnan(mod_traj1.X(:,j)),1);
    
    m_plot(mod_traj1.X(ids: ids+days_plt, j)-360, mod_traj1.Y(ids:ids+days_plt, j), 'linewidth',1, 'color',colors(2,:));

    m_plot(mod_traj1.X(ids, j)-360, mod_traj1.Y(ids, j), '.', 'Markersize',10,'color',colors(5,:));

end 

for j=1:4:240
    ids = find(~isnan(mod_traj2.X(:,j)),1);
    m_plot(mod_traj2.X(ids: ids+days_plt, j)-360, mod_traj2.Y(ids:ids+days_plt, j), 'linewidth',1, 'color',colors(1,:));

    m_plot(mod_traj2.X(ids, j)-360, mod_traj2.Y(ids, j), '.', 'Markersize',10,'color',colors(5,:));

end 
set(gca, 'fontsize', 16)


m_grid('linest', 'none')
saveas(gca,'../figures/intro_models.pdf')


