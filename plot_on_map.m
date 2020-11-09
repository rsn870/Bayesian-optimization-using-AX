%% Plot locations on a map

latlim_plot = [4 30];
lonlim_plot = [65 95];

m_proj('miller','lat',latlim_plot,'lon',lonlim_plot)
m_gshhs_f('save','gumby');

f1 = figure;
%set(gcf,'Position', [100, 100, 1200, 900]);
clf
set(gcf,'color','w')
hold on

m_proj('miller','lat',latlim_plot,'lon',lonlim_plot);

m_grid('box','fancy','fontsize',13)

% cmocean('algae');
%%cmocean('balance','zero');
m_usercoast('gumby','patch',[.7 .7 .7])

hold on

m_scatter(lonlist_k,latlist_k,20,betabias_k);
%clabel(C,H);
colorbar
%caxis([0 40])