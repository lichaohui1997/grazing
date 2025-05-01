
addpath(genpath('/Users/lichaohui/Desktop/calculation/grazingniche'));
addpath(genpath('/Users/lichaohui/Desktop/calculation/command'));

%% livestock density
load('grassland_env.mat')
% this livestockdensity.mat comes from the below processing of data
% load this into directory and you don't need to process the below n
% sections of code related to livestock. 
% this livestockdensity.mat data has gotten rid of the non-grassland
% livestocks,etc.
load('livestockdensity.mat')
% this AGB.mat data comes from the below section on grassland biomass.
load('AGB.mat')
load('aggregate_tas.mat')
load('niche.mat')% variable name: landuse_coupled1 (in run4)

load('aggregate_hurs.mat')
aggregate_hurs=aggregate_data;
load('aggregate_pr.mat')
aggregate_pr=aggregate_data;
load('aggregate_sfcWind.mat')
aggregate_sfcWind=aggregate_data;
load('aggregate_tas.mat')
aggregate_tas=aggregate_data;
%% modern livestock distribution niche
imagesc(aggregate_tas)
imagesc(resizecattle)
imagesc(aggregate_hurs)
colorbar

%%
nexttile;
hold on;
scatter(vector_pr, vector_cattle,'filled');
colormap(jet);
xlabel('Precipitation (mm)');
ylabel('Historical grazing intensity');
title('Niche Distribution - Precipitation');
defualtAxes()
% Define the x-coordinates for the grey area
xbars = [50 2627]; % x-coordinates for the grey area
yLimits = ylim; % Get the current y-axis limits
% Create the patch with EdgeColor set to 'none'
patch([xbars(1) xbars(1) xbars(2) xbars(2)], [yLimits(1) yLimits(2) yLimits(2) yLimits(1)], [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold off;

nexttile;
hold on;
scatter(vector_tas, vector_cattle,'filled');
colormap(jet);
xlabel('Mean Annual Temperature (^{o}C)');
ylabel('Historical grazing intensity');
title('Niche Distribution - Temperature');
defualtAxes()
% Define the x-coordinates for the grey area
xbars = [-3 29]; % x-coordinates for the grey area
yLimits = ylim; % Get the current y-axis limits
% Create the patch with EdgeColor set to 'none'
patch([xbars(1) xbars(1) xbars(2) xbars(2)], [yLimits(1) yLimits(2) yLimits(2) yLimits(1)], [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold off;

nexttile;
hold on;
scatter(vector_hurs, vector_cattle,'filled');
colormap(jet);
xlabel('Relative Humidity (%)');
ylabel('Historical grazing intensity');
title('Niche Distribution - Humidity');
defualtAxes()
% Define the x-coordinates for the grey area
xbars = [39 67]; % x-coordinates for the grey area
yLimits = ylim; % Get the current y-axis limits
% Create the patch with EdgeColor set to 'none'
patch([xbars(1) xbars(1) xbars(2) xbars(2)], [yLimits(1) yLimits(2) yLimits(2) yLimits(1)], [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold off;

nexttile;
hold on;
scatter(vector_sfcWind, vector_cattle,'filled');
colormap(jet);
xlabel('Near Surface Windspeed (m/s)');
ylabel('Historical grazing intensity');
title('Niche Distribution - Windspeed');
defualtAxes()
% Define the x-coordinates for the grey area
xbars = [1 6]; % x-coordinates for the grey area
yLimits = ylim; % Get the current y-axis limits
% Create the patch with EdgeColor set to 'none'
patch([xbars(1) xbars(1) xbars(2) xbars(2)], [yLimits(1) yLimits(2) yLimits(2) yLimits(1)], [0.8 0.8 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold off;

set(gcf, 'Position',  [478,192,778,687])
%%
imagesc(aggregate_sfcWind)
colorbar
caxis([0,10])
%% read population density data
[pop, Rpop] = readgeoraster('gpw_v4_population_density_rev11_2015_2pt5_min.tif')
Rpop = georefcells([-90,90],[-180,180],size(pop));
[resizedpop,resizedRpop] = georesize(pop,Rpop,1800/4320,"bilinear");
resizedpop(resizedpop<0)=0
imagesc(resizedpop)
colorbar
%% import cattle density and resize
% this section of code will take 30s to run
[cattle, Rcattle] = readgeoraster('/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/livestockDensity/cattle/Glb_Cattle_CC2006_AD.tif');

[allresizecattle,resizedRcattle] = georesize(cattle,Rcattle,1/12,"bilinear");
% all nan values are recorded as -1. this interferes with the calculation
allresizecattle(allresizecattle<0)=0

% grassland cattle
%resizecattle to equal to allresizecattle where grassland_env==0
resizecattle = allresizecattle; % First, make allresizecattle a copy of resizecattle
resizecattle(resizecattle>250) = 0;
resizecattle(resizedpop >20) = 0; % Then, set values to 0 where grassland_env is not 0
% cattlenum=sum(sum(resizecattle.*worldarea));
% cattlenumall=sum(sum(allresizecattle.*worldarea));


imagesc(allresizecattle)
colorbar

imagesc(resizecattle)
colorbar
%% R setting
R=georefcells([-90,90],[-180,180],size(resizecattle));
%% graph: cattle
figure1 = figure('WindowState','fullscreen');
sgtitle('Grassland cattle density', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1;addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar;
R=georefcells([-90,90],[-180,180],size(resizecattle));
geoshow(flipud(resizecattle), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/SI/cattle_density.png');

%% import sheep density and resize
[sheep, Rsheep] = readgeoraster('/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/livestockDensity/sheep/Glb_SHAD_2006.tif');
[allresizesheep,resizedRsheep] = georesize(sheep,Rsheep,1/12,"bilinear");
allresizesheep(allresizesheep<0) =0;

% keeping only the grassland sheep
resizesheep=allresizesheep;
resizesheep(resizedpop >20)=0
% sheepsnum=sum(sum(resizesheep.*worldarea));
% sheepsnumall=sum(sum(allresizesheep.*worldarea));
%% graph: sheep
figure1 = figure('WindowState','fullscreen');
sgtitle('Grassland sheep density', 'FontSize', 16);
custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1;addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
geoshow(flipud(resizesheep), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/SI/sheep_density.png');

%% import goats density and resize
[goats, Rgoats] = readgeoraster('/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/livestockDensity/goats/Glb_GTAD_2006.tif');
Rgoats=Rsheep;% somehow Rgoats is a empty array and does not contain georeferencing information
[allresizegoats,resizedRgoats] = georesize(goats,Rgoats,1/12,"bilinear");
allresizegoats(allresizegoats<0) =0;
allresizegoats(allresizegoats>250) = 250;
% keeping only the grassland goats
resizegoats=allresizegoats;
resizegoats(resizedpop >20)=0;
% goatsnum=sum(sum(resizegoats.*worldarea));
% goatsnumall=sum(sum(allresizegoats.*worldarea));

%% graph: goats
figure1 = figure('WindowState','fullscreen');
sgtitle('Grassland goats density', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1;addcolorplus(341)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
geoshow(flipud(resizegoats), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/SI/goats_density.png');
%% livestock unit
% Although I calculated this, I did not use this in my analysis
% I individually accounted for the feed of each of the animal types
% but If you want to use standardized animal unit you can use this
% equation. Note the parameters might be different. 
livestock=resizecattle+0.2*resizegoats+0.2*resizesheep
% %% livestock intensity graph
% figure1 = figure('WindowState','fullscreen');
% sgtitle('Grassland livestock density', 'FontSize', 16);
% %custom_colormap = [1 1 1; jet(255)]; 
% custom_colormap = [1 1 1;addcolorplus(341)];
% ax = worldmap('world')      
% setm(ax, 'FFaceColor', [1 1 1]);    
% %colormap(ax,custom_colormap);
% set(gcf, 'Colormap', custom_colormap);
% colorbar
% geoshow(flipud(livestock), R, 'DisplayType', 'texturemap');   
% load coastlines
% plotm(coastlat, coastlon, 'Color', 'black'); 
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/SI/livestock_density.png');
%% merged figure of livestock
% figure1 = figure('WindowState','fullscreen');
% sgtitle('Grassland Animal Density', 'FontSize', 16);

% Custom positions for each subplot

% % First map (Top-left)
% subplot(2,3,1);
% title('Grassland cattle density');
% worldmap('world');
% %setm(ax1, 'FFaceColor', [1 1 1]);
% custom_colormap = [1 1 1;addcolorplus(341)];
% set(gcf, 'Colormap', custom_colormap);
% geoshow(flipud(resizecattle), R, 'DisplayType', 'texturemap');
% load coastlines;
% plotm(coastlat, coastlon, 'Color', 'black');
% 
% % Second map (Top-middle)
% subplot(2,3,2);
% title('Grassland sheep density');
% worldmap('world');
% custom_colormap = [1 1 1;addcolorplus(341)];
% set(gcf, 'Colormap', custom_colormap);
% geoshow(flipud(resizesheep), R, 'DisplayType', 'texturemap');
% plotm(coastlat, coastlon, 'Color', 'black');
% 
% % Third map (Top-right)
% subplot(2,3,3);
% title('Grassland goats density');
% worldmap('world');
% custom_colormap = [1 1 1;addcolorplus(341)];
% set(gcf, 'Colormap', custom_colormap);
% geoshow(flipud(resizegoats), R, 'DisplayType', 'texturemap');
% plotm(coastlat, coastlon, 'Color', 'black');
% 
% % Fourth map (Bottom)
% subplot(2,3,[4,6]);
% title('Grassland livestock density');
% worldmap('world');
% custom_colormap = [1 1 1;addcolorplus(341)];
% set(gcf, 'Colormap', custom_colormap);
% colorbar;
% geoshow(flipud(livestock), R, 'DisplayType', 'texturemap');
% plotm(coastlat, coastlon, 'Color', 'black');
% 
% % Save the figure
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/manu/livestock_density_subplot.png');

% save the livestock data
save('./grazingniche/matdata/livestockdensity.mat','resizecattle','resizedRcattle','resizegoats','resizedRgoats','resizesheep','resizedRsheep','livestock')
load('livestockdensity.mat') %variables: 'resizecattle','resizedRcattle','resizegoats','resizedRgoats','resizesheep','resizedRsheep')

