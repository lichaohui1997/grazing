clear;clear all;
addpath(genpath('/Users/lichaohui/Desktop/calculation/grazingniche'));
addpath(genpath('/Users/lichaohui/Desktop/calculation/command'));

%% vector climate and grassland data
load('final_data.mat')
pr_final_vector=sum(pr_final,2)/7;
tas_final_vector=sum(tas_final,2)/7;
sfcWind_final_vector=sum(sfcWind_final,2)/7;
hurs_final_vector=sum(hurs_final,2)/7;
%% scatter plot e.g. tas-pr

subplot(2,3,1)
X=pr_final_vector;
Y=tas_final_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(min(X),max(X),n);
YList=linspace(min(Y),max(Y),n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);

% 
scatter(X,Y,'filled','CData',H);
xlabel('Precipitation (mm)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Niche Distribution'; 'Precipitation and Temperature'});
defualtAxes()

subplot(2,3,2)
X=hurs_final_vector;
Y=tas_final_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(min(X),max(X),n);
YList=linspace(min(Y),max(Y),n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);

% 
scatter(X,Y,'filled','CData',H);
xlabel('Relative humidity (%)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Niche Distribution'; 'Humidity and Temperature'});
defualtAxes()

subplot(2,3,3)
X=sfcWind_final_vector;
Y=tas_final_vector;
% 横竖分割一百格计算核密度
n=100;
XList=linspace(min(X),max(X),n);
YList=linspace(min(Y),max(Y),n);
[XMesh,YMesh]=meshgrid(XList,YList);
F=ksdensity([X,Y],[XMesh(:),YMesh(:)]);
ZMesh=reshape(F,size(XMesh));
H=interp2(XMesh,YMesh,ZMesh,X,Y);

% 
scatter(X,Y,'filled','CData',H);
xlabel('Near Surface Windspeed (m/s)');
ylabel('Mean Annual Temperature (^{o}C)');
title({'Niche Distribution'; 'Windspeed and Temperature'});
defualtAxes()

% Save the figure as a SVG and specify the size
set(gcf, 'Position',  [751,163,1092,753])

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/subplot_niche_twowayscatter.svg');

%% aggregate climate data (matrix form) from run4

load('aggregate_hurs.mat')
aggregate_hurs=aggregate_data;
load('aggregate_pr.mat')
aggregate_pr=aggregate_data;
load('aggregate_sfcWind.mat')
aggregate_sfcWind=aggregate_data;
load('aggregate_tas.mat')
aggregate_tas=aggregate_data;

%% grassland data
load('grassland_env.mat')
%% Niche

cond1_pr = (aggregate_pr >= 50 & aggregate_pr <= 2627);
cond1_tas = (aggregate_tas >= -3 & aggregate_tas <= 29);
cond1_sfcWind = (aggregate_sfcWind >= 1 & aggregate_sfcWind <= 6);
cond_hurs = (aggregate_hurs >= 39 & aggregate_hurs <= 67);
cond1_landuse = (grassland_env>0);

% Combine the conditions to create the niche map
niche1 = cond1_pr & cond1_tas & cond1_sfcWind & cond_hurs & cond1_landuse;

% Create the coupled landuse map
landuse_coupled1 = double(niche1) .* grassland_env;

% producing the graph of grazing niche
R = georefcells([-90,90],[-180,180],size(landuse_coupled1));

figure1 = figure;
%sgtitle('Niche map', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
%custom_colormap = [1 1 1; addcolorplus(332)];
custom_colormap = [1 1 1; addcolorplus(309)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
c = colorbar;  
c.Ruler.TickLabelFormat='%g%%';
geoshow(flipud(landuse_coupled1), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
%set(gcf, 'Position',  [584,449,684,537])
set(gcf, 'Position',  [361,614,899,295])


saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/niche_general.svg');
save('./grazingniche/matdata/niche.mat',"landuse_coupled1")

%% grassland figure
% to run this you need to add the command file into directory
figure1 = figure;
sgtitle('Grassland map', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(332)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
c.Ruler.TickLabelFormat='%g%%';
R = georefcells([-90,90],[0,360],size(grassland_env));
geoshow(flipud(grassland_env), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
set(gcf, 'Position',  [584,449,684,537])

% set(gcf,'renderer','painters');
% it seems matlab is unable to produce real svg with this map
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/grassland.svg');
