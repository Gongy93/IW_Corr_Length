clc
clear
close all

%% Plot the correlation length scales for internal waves 
% (see Mazloff et al. 2018, Gong et al. 2021)

% Aims:
% Plot the corrlation length scales for internal waves, 
% including: L_min: semiminor axis (corresponding to crest widths)
%            L_maj: semiminor axis (corresponding to crest lengths)
%            sida: ellipse orientation (corresponding to propagation directions)

% Input:
% Matlab file: "Corr_ellipse_params.mat".
% dx: horizontal resolution [unit: m]
% x0: selected box size, e.g. 20 km * 20 km
% x_interval = 5: do the loop every 5 grids (= 0.5 km)

% Output:
% The Image of internal wave properties in space.

% Note:
% Before plot the L_min and L_maj, we need to make sure "L_maj > L_min" 
% at each point, otherwise swap their value.

%  written by Gong Yankun in February 2021

%% Load the correlation ellipse parameters
load('Corr_ellipse_params.mat');

%% Swap semimajor axis and semiminor axis if it is required

% Params setup (same as those in the file "get1_Var_correlation.m")
dx = 100; % horizontal resolution [unit: m]
x0 = 20000; % selected box size: 20 km * 20 km
x_interval = 5; % do the loop every 5 grids (= 0.5 km)

% Define the starting point and ending point to do the loop
start_lon = x0/dx; 
end_lon = size(a1,1) - x0/dx;
start_lat = x0/dx;
end_lat = size(a1,2) - x0/dx;

MIN_a = min_a;
MAX_a = max_a;

for lo = start_lon :x_interval: end_lon
   for la = start_lat :x_interval: end_lat

   if MIN_a(lo,la) > MAX_a(lo,la)
       min_a(lo,la) = MAX_a(lo,la);
       max_a(lo,la) = MIN_a(lo,la);

 % As sida is the orientation of the ellipse, which means the orientation of 
 % semimajor axis. However, the internal waves propagate in the direction
 % of semiminor axis. We therefore rotate it 90 degree.
 
   if sida(lo,la)>0
       sida(lo,la) = sida(lo,la)-90;
   else if sida(lo,la)<0
       sida(lo,la) = sida(lo,la)+90;   
   end
       
   end
   
   end
end
end

rate_a = min_a ./ max_a;

%% Load the model grid information and topography dataset for plotting
Xgrid = load('Xgrid_meters.txt')/1000; % unit: km
Ygrid = load('Ygrid_meters.txt')/1000;

prec='real*8';
ieee='b';
fid=fopen('RowleyShoals.topo.3d','r',ieee);
depth=fread(fid,prec);
fclose(fid);
depth = reshape(depth,size(a1,1),size(a1,2));

%% Here are two figures below, first for the correlation ellipses, 
% the second for internal wave properties (L_min, L_maj, sida)
colormap_ZeorAtBeginV2 = dlmread('colormap_ZeorAtBegin.txt'); 

%% Figure 1:
% Note that the factor of 0.4723 is for transfering L_maj to crest length
% and L_min to crest width (see Appendix in Gong et al. 2021)
max_a(max_a*0.4723/10>30) = 299/0.4723;
min_a(min_a*0.4723/10>30) = 299/0.4723;

figure('visible','off','position',[100 100 1000 1600]);
% a1
a1 = axes('position',[0.12 0.1 0.8 0.82]);
for lo = start_lon:40:end_lon
    for la = start_lat:40:end_lat  
       hold on; Ellipse_new(max_a(lo,la)/250,min_a(lo,la)/250,...
           sida(lo,la)*2*pi/360,Xgrid(lo),Ygrid(la),'b'); % 
% Ellipse_new is a toolbox for ellipse plotting, written by D.G. Long

    end
end

hold on; contour(Xgrid,Ygrid,depth',...
  'LevelList',[0],'linecolor',[0 0 0],'linewidth',2)  
hold on
[C h] = contour(Xgrid,Ygrid,-depth',...
  'LevelList',-[-50 -100 -150 -200 -250 -300 -400 -500 -600 -800 -1000 -1200 -1600],'linecolor',[0.3 0.3 0.3],'linewidth',1.8); 
clabel(C,h,'labelspacing',2000,'fontsize',13,'color','k');
box on

xlabel('x [km]','fontsize',15);
ylabel('y [km]','fontsize',15);
title('Correlation Ellipses of Internal Waves','fontsize',18);
xlim([Xgrid(start_lon/2) Xgrid(end-start_lon/2)]);
ylim([Ygrid(start_lat/2) Ygrid(end-start_lat/2)]); 
box on

set(a1,'fontsize',15,'linewidth',2);
set(gcf,'paperpositionmode','auto');
saveas(gcf,'Corr_ellipse_IWs.png','png');
close  

%% Figure 2:
% 4 panels plot: Patchy colorbar

% Mask out the shallow water zone (<200 m)
max_a(depth>-200) = nan;
min_a(depth>-200) = nan;
rate_a(depth>-200) = nan;
sida(depth>-200) = nan;

figure('visible','off','position',[100 100 780 750]);

% a1 : Major axis 
a1 = axes('position',[0.08 0.55 0.38 0.4]);
pcolor(Xgrid(start_lon:x_interval:end_lon),...
    Ygrid(start_lat:x_interval:1200),...
    max_a(start_lon:x_interval:end_lon,start_lat:x_interval:1200)'*0.4723/10);shading interp
caxis([.0 30]); 
colormap(a1,colormap_ZeorAtBeginV2(5:5:end-40,:));%colorbar
cbarid = colorbar('position',[0.5-0.02 0.55 0.01 0.4]);
set(get(cbarid,'Title'),'string',' km ','fontsize',12); 
hold on; contour(Xgrid,Ygrid,depth',...
  'LevelList',[0],'linecolor',[0 0 0],'linewidth',1.8)  
hold on
[C h] = contour(Xgrid(100:1400),Ygrid(100:1400),-depth(100:1400,100:1400)',20,...
        'LevelList',[150 200 250 300 350 400 500 600 800 1000 1600 2000 2500],...
      'linecolor',[0.5 0.5 0.5],'linewidth',1.3);
clabel(C,h,'labelspacing',2000,'fontsize',9,'color',[.5 .5 .5]); 
box on
title(' L_{maj} ','fontsize',12);
ylabel('y [km]','fontsize',12);
xlim([Xgrid(150) Xgrid(end-150)]);
ylim([Ygrid(150) Ygrid(end-750)]); 

% a2 : Minor axis 
a2 = axes('position',[0.54 0.55 0.38 0.4]);
pcolor(Xgrid(start_lon:x_interval:end_lon),...
    Ygrid(start_lat:x_interval:1200),...
    min_a(start_lon:x_interval:end_lon,start_lat:x_interval:1200)'*0.4723/10);shading interp
caxis([.0 30]); 
colormap(a2,colormap_ZeorAtBeginV2(5:5:end-40,:));
cbarid = colorbar('position',[0.94 0.55 0.01 0.4]);
set(get(cbarid,'Title'),'string',' km ','fontsize',12); 
hold on; contour(Xgrid,Ygrid,depth',...
  'LevelList',[0],'linecolor',[0 0 0],'linewidth',1.8)  
hold on
[C h] = contour(Xgrid(100:1400),Ygrid(100:1400),-depth(100:1400,100:1400)',20,...
        'LevelList',[150 200 250 300 350 400 500 600 800 1000 1600 2000 2500],...
      'linecolor',[0.5 0.5 0.5],'linewidth',1.3);
 clabel(C,h,'labelspacing',2000,'fontsize',9,'color',[.5 .5 .5]);   
box on
title(' L_{min} ','fontsize',12);
xlim([Xgrid(150) Xgrid(end-150)]);
ylim([Ygrid(150) Ygrid(end-750)]); 

% a3 : Ratio of min_a to max_a
rate_a(rate_a>1) = 0.99;
a3 = axes('position',[0.08 0.08 0.38 0.4]);
pcolor(Xgrid(start_lon:x_interval:end_lon),...
    Ygrid(start_lat:x_interval:1200),...
    rate_a(start_lon:x_interval:end_lon,start_lat:x_interval:1200)');shading interp
caxis([0 1]); 
colormap(a3,parula(30));
cbarid = colorbar('position',[0.5-0.02 0.08 0.01 0.4]);
set(get(cbarid,'Title'),'string',' km km^-^1 ','fontsize',12); 
hold on; contour(Xgrid,Ygrid,depth',...
  'LevelList',[0],'linecolor',[0 0 0],'linewidth',1.8)  
hold on
[C h] = contour(Xgrid(100:1400),Ygrid(100:1400),-depth(100:1400,100:1400)',20,...
        'LevelList',[150 200 250 300 350 400 500 600 800 1000 1600 2000 2500],...
      'linecolor',[0.5 0.5 0.5],'linewidth',1.3);
  clabel(C,h,'labelspacing',2000,'fontsize',9,'color',[.5 .5 .5]);
box on
title(' L_{min} / L_{maj} ','fontsize',12);
xlabel('x [km]','fontsize',12);
ylabel('y [km]','fontsize',12);
xlim([Xgrid(150) Xgrid(end-150)]);
ylim([Ygrid(150) Ygrid(end-750)]); 

% a4 : Orientation 
% plot the vectors to show the propagation directions of IWs
nor_y0 = - ones(length(Xgrid),length(Ygrid));
nor_x0 = - nor_y0 .* tand(sida);
nor_x = nor_x0 ./ sqrt(nor_x0.^2 + nor_y0.^2);
nor_y = nor_y0 ./ sqrt(nor_x0.^2 + nor_y0.^2); clear nor_x0 nor_y0
[mesh_x,mesh_y] = meshgrid(Xgrid,Ygrid);
mesh_x = mesh_x'; mesh_y = mesh_y';

% Exclude the shallow water region
sida(depth>-200) = nan;
nor_x(depth>-200) = nan;
nor_y(depth>-200) = nan;

a4 = axes('position',[0.54 0.08 0.38 0.4]);
pcolor(Xgrid(start_lon:x_interval:end_lon),...
    Ygrid(start_lat:x_interval:1200),...
    sida(start_lon:x_interval:end_lon,start_lat:x_interval:1200)');shading interp
caxis([-90 90]); 
colormap(a4,flipud(othercolor('RdYlBu6',30)));
cbarid = colorbar('position',[0.94 0.08 0.01 0.4]);
set(get(cbarid,'Title'),'string',' deg ','fontsize',12); 
hold on; contour(Xgrid,Ygrid,depth',...
  'LevelList',[0],'linecolor',[0 0 0],'linewidth',1.8)  
hold on
[C h] = contour(Xgrid,Ygrid,-depth',...
  'LevelList',-[-50 -100 -150 -200 -250 -300 -350 -400 -500 -600 -800 -1000 -1200 -1600],...
  'linecolor',[0.5 0.5 0.5],'linewidth',1.2);
hold on;
 quiver(mesh_x(220:50:end,220:50:1200),mesh_y(220:50:end,220:50:1200),...
   nor_x(220:50:end,220:50:1200),nor_y(220:50:end,220:50:1200),...
  0.8,'Color',[0 0 0],'linewidth',1) ;   
box on
title(' Orientation \theta - 90^o','fontsize',12);
xlabel('x [km]','fontsize',12);
xlim([Xgrid(150) Xgrid(end-150)]);
ylim([Ygrid(150) Ygrid(end-750)]); 

set(a1,'fontsize',11,'linewidth',2,'xticklabel',[]);
set(a2,'fontsize',11,'linewidth',2,'xticklabel',[],'yticklabel',[]);
set(a3,'fontsize',11,'linewidth',2);
set(a4,'fontsize',11,'linewidth',2,'yticklabel',[]);

set(gcf,'paperpositionmode','auto');
saveas(gcf,'Corr_lengthscale_IWs.png','png');
close  
