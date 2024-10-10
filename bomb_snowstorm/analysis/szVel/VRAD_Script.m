clc; clear all; close all;

indir='/scr/cirrus1/rsfdata/projects/nexrad/matFiles/';
figdir='/scr/cirrus1/rsfdata/projects/nexrad/figures/szComp/';

load([indir,'VRAD_KFTG_Case.mat']);


axv = [-150 210 -180 180];


figure; set(gcf,'units','normalized','outerposition',[0.12 0.08 0.50 0.88]);


subplot(2,2,1);
pcolor(xx,yy,z.*thr); 
shading flat; axis equal; caxis([-20 80]); 
xlabel('Distance from Radar E/W (km)'); 
ylabel('Distance from Radar N/S (km)');
title('Reflectivity (dBZ)'); axis(axv); 
bmapover(gca,[1 1 1 zeros(1,9)],radar_name);
h1 = colormap(gca,czmap); 
hb = colorbar; set(gcf,'color','w'); 
set(get(hb,'ylabel'),'string','dBZ');
set(gca,'linewidth',2,'Fontsize',16,'FontName','Times New Roman');
hold on; plot([0 0], [0 0], 'k.','markersize',20);


subplot(2,2,2);
pcolor(xx,yy,v1.*thr); 
shading flat; axis equal; caxis([-30 30]); 
h2 = colormap(gca,rgmap); cb1 = colorbar;
xlabel('Distance from Radar E/W (km)'); 
ylabel('Distance from Radar N/S (km)');
title('CD Velocities'); axis(axv); 
bmapover(gca,[1 1 1 zeros(1,9)],radar_name);
set(get(cb1,'ylabel'),'string','m s^{-1}');
set(gca,'linewidth',2,'Fontsize',16,'FontName','Times New Roman');
hold on; plot([0 0], [0 0], 'k.','markersize',20);


subplot(2,2,3);
pcolor(xx,yy,v2.*thr); 
shading flat; axis equal; caxis([-9 9]); 
colormap(gca,rgmap); cb = colorbar;
xlabel('Distance from Radar E/W (km)'); 
ylabel('Distance from Radar N/S (km)');
bmapover(gca,[1 1 1 zeros(1,9)],radar_name);
title('CS Velocities'); axis(axv);
set(get(cb,'ylabel'),'string','m s^{-1}');
set(gca,'linewidth',2,'Fontsize',16,'FontName','Times New Roman');
hold on; plot([0 0], [0 0], 'k.','markersize',20);


subplot(2,2,4);
pcolor(xx,yy,vrad.*thr); 
shading flat; axis equal; caxis([-30 30]); 
h3 = colormap(gca,rgmap); cb2 = colorbar;
xlabel('Distance from Radar E/W (km)'); 
ylabel('Distance from Radar N/S (km)');
title('VRAD Velocities'); axis(axv); 
set(get(cb2,'ylabel'),'string','m s^{-1}');
set(gcf,'color','w'); bmapover(gca,[1 1 1 zeros(1,9)],radar_name);
set(gca,'linewidth',2,'Fontsize',16,'FontName','Times New Roman');
hold on; plot([0 0], [0 0], 'k.','markersize',20);


ax = findobj('Type','axes','Tag','');
ax_index = length(ax);
if ax_index > 1
    linkaxes(flipud(ax),'xy');
end

print([figdir,'VEL_VRAD.png'],'-dpng','-r0');