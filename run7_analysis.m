addpath(genpath('/Users/lichaohui/Desktop/calculation/grazingniche'));
addpath(genpath('/Users/lichaohui/Desktop/calculation/command'));

%% import data
% in order for this .m file to work, you would need to import:
load('grassland_env.mat')
load('livestockdensity.mat') % comes with these variables: 'resizecattle','resizedRcattle','resizegoats','resizedRgoats','resizesheep','resizedRsheep'
load('niche.mat') %variable name: landuse_coupled1 (in run4)
load('futureniche.mat')
load('futureland.mat')
load('livestock.mat')


%% continental analysis
x=ncinfo('/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc')
lon_nc = ncread('/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc', 'lon');
lat_nc = ncread('/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc', 'lat');
[lon_nc, lat_nc] = meshgrid(lon_nc, lat_nc);

data_world = ncread('/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc', 'm_world');
data_AUS = ncread('/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc', 'm_AUS');
data_CHN = ncread('/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc', 'm_CHN');

Rmask=georefcells([-90,90],[-180,180],size(data_world));

[resizeworld,resizedRworld] = georesize(data_world,Rmask,5,"bilinear");
[resizeAUS,resizedRAUS] = georesize(data_AUS,Rmask,5,"bilinear");
[resizeCHN,resizedRworld] = georesize(data_CHN,Rmask,5,"bilinear");

% the new shortened code
% Define path to your file
file_path = '/Users/lichaohui/Desktop/calculation/command/countrymasks-master/countrymasks.nc';

% countries = {'m_AFG', 'm_ALB', 'm_DZA', 'm_AND', 'm_AGO', 'm_ATG', 'm_ARG', 'm_ARM', ...
%              'm_AUS', 'm_AUT', 'm_AZE', 'm_BHS', 'm_BHR', 'm_BGD', 'm_BRB', 'm_BLR', ...
%              'm_BEL', 'm_BLZ', 'm_BEN', 'm_BTN', 'm_BOL', 'm_BIH', 'm_BWA', 'm_BRA', ...
%              'm_BRN', 'm_BGR', 'm_BFA', 'm_BDI', 'm_KHM', 'm_CMR', 'm_CAN', 'm_CPV', ...
%              'm_CSID', 'm_CYM', 'm_CAF', 'm_TCD', 'm_CHL', 'm_CHN', 'm_COL', 'm_COM', ...
%              'm_COG', 'm_CRI', 'm_HRV', 'm_CUB', 'm_CYP', 'm_CZE', 'm_CIV', 'm_PRK', ...
%              'm_COD', 'm_DNK', 'm_DJI', 'm_DMA', 'm_DOM', 'm_ECU', 'm_EGY', 'm_SLV', ...
%              'm_GNQ', 'm_ERI', 'm_EST', 'm_ETH', 'm_FLK', 'm_FRO', 'm_FJI', 'm_FIN', ...
%              'm_FRA', 'm_GUF', 'm_PYF', 'm_ATF', 'm_GAB', 'm_GMB', 'm_GEO', 'm_DEU', ...
%              'm_GHA', 'm_GRC', 'm_GRL', 'm_GRD', 'm_GLP', 'm_GUM', 'm_GTM', 'm_GIN', ...
%              'm_GNB', 'm_GUY', 'm_HTI', 'm_HMD', 'm_HND', 'm_HKG', 'm_HUN', 'm_ISL', ...
%              'm_IND', 'm_IOSID', 'm_IDN', 'm_IRN', 'm_IRQ', 'm_IRL', 'm_IMN', 'm_ISR', ...
%              'm_ITA', 'm_JAM', 'm_JKX', 'm_JPN', 'm_JOR', 'm_KAZ', 'm_KEN', 'm_KIR', ...
%              'm_KWT', 'm_KGZ', 'm_LAO', 'm_LVA', 'm_LBN', 'm_LSO', 'm_LBR', 'm_LBY', ...
%              'm_LTU', 'm_LUX', 'm_MDG', 'm_MWI', 'm_MYS', 'm_MLI', 'm_MLT', 'm_MTQ', ...
%              'm_MRT', 'm_MUS', 'm_MYT', 'm_MEX', 'm_FSM', 'm_MDA', 'm_MNG', 'm_MNE', ...
%              'm_MAR', 'm_MOZ', 'm_MMR', 'm_NAM', 'm_NPL', 'm_NLD', 'm_ANT', 'm_NCL', ...
%              'm_NZL', 'm_NIC', 'm_NER', 'm_NGA', 'm_NIU', 'm_NOR', 'm_OMN', 'm_PSID', ...
%              'm_PAK', 'm_PLW', 'm_PSE', 'm_PAN', 'm_PNG', 'm_PRY', 'm_PER', 'm_PHL', ...
%              'm_POL', 'm_PRT', 'm_PRI', 'm_QAT', 'm_KOR', 'm_ROU', 'm_RUS', 'm_RWA', ...
%              'm_REU', 'm_LCA', 'm_SPM', 'm_VCT', 'm_WSM', 'm_STP', 'm_SAU', 'm_SEN', ...
%              'm_SRB', 'm_SLE', 'm_SGP', 'm_SVK', 'm_SVN', 'm_SLB', 'm_SOM', 'm_ZAF', ...
%              'm_SGS', 'm_SSD', 'm_ESP', 'm_LKA', 'm_SDN', 'm_SUR', 'm_SJM', 'm_SWZ', ...
%              'm_SWE', 'm_CHE', 'm_SYR', 'm_TWN', 'm_TJK', 'm_THA', 'm_MKD', 'm_TLS', ...
%              'm_TGO', 'm_TON', 'm_TTO', 'm_TUN', 'm_TUR', 'm_TKM', 'm_GBR', 'm_UGA', ...
%              'm_UKR', 'm_ARE', 'm_TZA', 'm_VIR', 'm_USA', 'm_URY', 'm_UZB', 'm_VUT', ...
%              'm_VEN', 'm_VNM', 'm_ESH', 'm_YEM', 'm_ZMB', 'm_ZWE', 'm_world'};

countries = {x.Variables(3:217).Name};

% Create a structure to hold all country data
country_data = struct();

% Read lat and lon, and create a meshgrid
lon = ncread(file_path, 'lon');
lat = ncread(file_path, 'lat');
[lon, lat] = meshgrid(lon, lat);
Rmask = georefcells([-90,90],[-180,180], size(ncread(file_path, 'm_world')));

% Loop through countries, read and resize each mask
for i = 1:length(countries)
    country_name = countries{i};
    
    % Read country data
    country_data_raw = ncread(file_path, country_name);
    
    % Resize the country data
    [resized_data, ~] = georesize(country_data_raw, Rmask, 5, "bilinear");
    
    % Store in the structure
    country_data.(country_name) = resized_data';
end

Asia={'m_CHN','m_JPN','m_KOR'}
Asiamask=country_data.m_CHN+country_data.m_JPN+country_data.m_KOR
imagesc(Asiamask)
colorbar

% mask for continents
% asia
Asia = {'m_AFG', 'm_ARM', 'm_AZE', 'm_BHR', 'm_BGD', 'm_BTN', 'm_BRN', 'm_KHM', 'm_CHN', 'm_CYM',... 
        'm_GEO', 'm_IND', 'm_IDN', 'm_IRN', 'm_IRQ', 'm_ISR', 'm_JPN', 'm_JOR', 'm_KAZ', 'm_KWT', ...
        'm_KGZ', 'm_LAO', 'm_LBN', 'm_MYS', 'm_MNG', 'm_MMR', 'm_NPL', 'm_PRK', 'm_OMN', ...
        'm_PAK', 'm_PHL', 'm_QAT', 'm_RUS', 'm_SAU', 'm_SGP', 'm_KOR', 'm_LKA', 'm_SYR', 'm_TWN', ...
        'm_TJK', 'm_THA', 'm_TUR', 'm_TKM', 'm_ARE', 'm_UZB', 'm_VNM', 'm_YEM'};

Asiamask = country_data.m_AFG;
for i = 2:length(Asia)
    Asiamask = Asiamask + country_data.(Asia{i});
end

imagesc(Asiamask)
colorbar
% for europe
Europe = {'m_ALB', 'm_AND', 'm_AUT', 'm_BEL', 'm_BIH', 'm_BGR', 'm_HRV', 'm_CYP', 'm_CZE', 'm_DNK',... 
          'm_EST', 'm_FIN', 'm_FRA', 'm_DEU', 'm_GRC', 'm_HUN', 'm_ISL', 'm_IRL', 'm_ITA', ...
          'm_LVA', 'm_LTU', 'm_LUX', 'm_MLT', 'm_MDA', 'm_MNE', 'm_NLD', 'm_MKD', ...
          'm_NOR', 'm_POL', 'm_PRT', 'm_ROU', 'm_SRB', 'm_SVK', 'm_SVN', 'm_ESP', ...
          'm_SWE', 'm_CHE', 'm_UKR', 'm_GBR'};

Europemask = country_data.m_ALB;
for i = 2:length(Europe)
    Europemask = Europemask + country_data.(Europe{i});
end

imagesc(Europemask)
colorbar

% now for south america
SouthAmerica = {'m_ARG', 'm_BOL', 'm_BRA', 'm_CHL', 'm_COL', 'm_ECU', 'm_GUY', 'm_PRY', 'm_PER', 'm_SUR', 'm_URY', 'm_VEN'};
SouthAmericamask = country_data.m_ARG;
for i = 2:length(SouthAmerica)
    SouthAmericamask = SouthAmericamask + country_data.(SouthAmerica{i});
end
imagesc(SouthAmericamask)
colorbar

% north america
NorthAmerica = {'m_ATG', 'm_BHS', 'm_BRB', 'm_BEL', 'm_CAN', 'm_CYM', 'm_CRI', 'm_CUB', 'm_DMA', 'm_DOM', 'm_SLV', 'm_GRL', 'm_GRD', 'm_GLP', 'm_GTM', 'm_HND', 'm_JAM', 'm_MEX', 'm_MTQ', 'm_NIC', 'm_PAN', 'm_PRI', 'm_SPM', 'm_LCA', 'm_VCT', 'm_TTO', 'm_USA', 'm_VIR'};
NorthAmericamask = country_data.m_ATG;
for i = 2:length(NorthAmerica)
    NorthAmericamask = NorthAmericamask + country_data.(NorthAmerica{i});
end
imagesc(NorthAmericamask)
colorbar
% africa
Africa = {'m_DZA', 'm_COD','m_SDN', 'm_AGO', 'm_BEN', 'm_BWA', 'm_BFA', 'm_BDI', 'm_CMR', 'm_CPV', 'm_CAF', 'm_TCD', 'm_COM', 'm_COG', 'm_CIV', 'm_DJI', 'm_EGY', 'm_GNQ', 'm_ERI', 'm_ETH', 'm_GAB', 'm_GMB', 'm_GHA', 'm_GIN', 'm_GNB', 'm_KEN', 'm_LSO', 'm_LBR', 'm_LBY', 'm_MDG', 'm_MWI', 'm_MLI', 'm_MRT', 'm_MUS', 'm_MAR', 'm_MOZ', 'm_NAM', 'm_NER', 'm_NGA', 'm_REU', 'm_RWA', 'm_STP', 'm_SEN', 'm_SLE', 'm_SOM', 'm_ZAF', 'm_SSD', 'm_ESH', 'm_SWZ', 'm_TZA', 'm_TGO', 'm_TUN', 'm_UGA', 'm_ZMB', 'm_ZWE'};
Africamask = country_data.m_DZA;
for i = 2:length(Africa)
    Africamask = Africamask + country_data.(Africa{i});
end
imagesc(Africamask)
colorbar

% Oceania
Oceania = {'m_AUS', 'm_FJI', 'm_KIR', 'm_FSM', 'm_NCL', 'm_NZL', 'm_NIU', 'm_PLW', 'm_PNG', 'm_WSM', 'm_SLB', 'm_TON', 'm_VUT'};

Oceaniamask = country_data.m_AUS;
for i = 2:length(Oceania)
    Oceaniamask = Oceaniamask + country_data.(Oceania{i});
end
imagesc(Oceaniamask)
colorbar

% convert the masks to double so that it's easier to manipulate them later
% then not only can they be used as logical layers we can directly multiply
% them with matricies, so we don't have to resort to too many conditions
% because the 0 is able to cancle any variable we want

Asiamask=double(Asiamask);
Europemask=double(Europemask);
Oceaniamask=double(Oceaniamask);
Africamask=double(Africamask);
NorthAmericamask=double(NorthAmericamask);
SouthAmericamask=double(SouthAmericamask);

save('./grazingniche/matdata/continentmask.mat','Asiamask','Europemask','Oceaniamask','Africamask','NorthAmericamask','SouthAmericamask')

%% calculating area

R = 6371; % Earth's radius in km
delta_lambda = 0.1 * pi / 180; % change in longitude in radians
phi1 = -pi/2; % starting from -90 degrees
phi2 = phi1 + 0.1 * pi / 180; % add 0.1 degree in radians

area = R^2 * delta_lambda * (sin(phi2) - sin(phi1));
disp(['The area of the (1,1) grid cell is approximately ', num2str(area), ' km^2']);

R = 6371; % Earth's radius in km
delta_lambda = 0.1 * pi / 180; % change in longitude in radians

% Initialize a vector to store the areas for each latitude.
areas = zeros(1800, 1);

% Loop over each latitude step
for i = 1:1800
    phi1 = (-pi/2) + (i-1) * 0.1 * pi / 180; % starting latitude for this step
    phi2 = phi1 + 0.1 * pi / 180; % ending latitude for this step

    areas(i) = R^2 * delta_lambda * (sin(phi2) - sin(phi1));
end

disp(areas);

worldarea=repmat(areas,[1,3600])

save('./grazingniche/matdata/worldarea.mat','worldarea');     
%%
earthsurface=sum(sum(worldarea))
% validated. very close to the googled value: 509 600 000 square km
%% land surface area
landsurfacearea=sum(sum(worldarea(resizeworld==1)));
Asiaarea=sum(sum(worldarea(Asiamask==1)));
Africaarea=sum(sum(worldarea(Africamask==1)));
SouthAmericAarea=sum(sum(worldarea(SouthAmericamask==1)));
NorthAmericaArea=sum(sum(worldarea(NorthAmericamask==1)));
Oceaniaarea=sum(sum(worldarea(Oceaniamask==1)));
Europearea=sum(sum(worldarea(Europemask==1)));
Landarea=[Asiaarea;Africaarea;SouthAmericAarea;NorthAmericaArea;Oceaniaarea;Europearea];
%% grassland area
% this section of code take 15s to run
% how much grassland are there in the world?
grasslandarea=grassland_env.*worldarea./100
grasslandarea_world=sum(sum(grasslandarea))

%how much grassland are there in each continent?
grasslandarea_Asia=sum(sum(grasslandarea(Asiamask==1)))
grasslandarea_Africa=sum(sum(grasslandarea(Africamask==1)))
grasslandarea_Oceania=sum(sum(grasslandarea(Oceaniamask==1)))
grasslandarea_SouthAmerica=sum(sum(grasslandarea(SouthAmericamask==1)))
grasslandarea_NorthAmerica=sum(sum(grasslandarea(NorthAmericamask==1)))
grasslandarea_Europe=sum(sum(grasslandarea(Europemask==1)))

grasslandarea=[grasslandarea_Asia;grasslandarea_Africa;grasslandarea_SouthAmerica;grasslandarea_NorthAmerica;grasslandarea_Oceania;grasslandarea_Europe]

% how much percent does that take?
grasslandprc=grasslandarea_world/landsurfacearea;
grasslandprc_Asia=grasslandarea_Asia/Asiaarea;
grasslandprc_Oceania=grasslandarea_Oceania/Oceaniaarea;
grasslandprc_NorthAmerica=grasslandarea_NorthAmerica/NorthAmericaArea;
grasslandprc_SouthAmerica=grasslandarea_SouthAmerica/SouthAmericAarea;
grasslandprc_Africa=grasslandarea_Africa/Africaarea;
grasslandprc_Europe=grasslandarea_Europe/Europearea;

%% niche area
nicheland_world=landuse_coupled1.*worldarea./100;
nichelandarea_world=sum(sum(nicheland_world));

nichelandarea_world_prct=nichelandarea_world./landsurfacearea

nichelandarea_Asia=sum(sum(nicheland_world(Asiamask==1)));
nichelandarea_Africa=sum(sum(nicheland_world(Africamask==1)));
nichelandarea_SouthAmerica=sum(sum(nicheland_world(SouthAmericamask==1)));
nichelandarea_NorthAmerica=sum(sum(nicheland_world(NorthAmericamask==1)));
nichelandarea_Oceania=sum(sum(nicheland_world(Oceaniamask==1)));
nichelandarea_Europe=sum(sum(nicheland_world(Europemask==1)));

save('./grazingniche/matdata/nichelandarea.mat','nichelandarea_Asia','nichelandarea_Africa','nichelandarea_SouthAmerica','nichelandarea_NorthAmerica','nichelandarea_Oceania','nichelandarea_Europe');     


nichelandarea=[nichelandarea_Asia;nichelandarea_Africa;nichelandarea_SouthAmerica;nichelandarea_NorthAmerica;nichelandarea_Oceania;nichelandarea_Europe]
nonnichelandarea=grasslandarea-nichelandarea

figure_nichelandarea_value=[nichelandarea,nonnichelandarea];

save('./grazingniche/matdata/figure_nichelandarea_value.mat','figure_nichelandarea_value');     

%% bar figure of niche area 
figure_nichelandarea_caption={'Asia';'Africa';'South America';'North America';'Oceania';'Europe'}

bar(figure_nichelandarea_value);
% Set x-axis labels
set(gca, 'xticklabel', figure_nichelandarea_caption);
% Add a title and labels
title('GN grassland vs Non-GN grassland');
xlabel('Region');
ylabel('Grassland area');
legend('GN area', 'Non-GN area');%% how much of the unusable grassland are overgrazed?

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/nichelandareabar.svg');

nicheprc_Asia=nichelandarea_Asia/grasslandarea_Asia;
nicheprc_Africa=nichelandarea_Africa/grasslandarea_Africa;
nicheprc_SouthAmerica=nichelandarea_SouthAmerica/grasslandarea_SouthAmerica;
nicheprc_NorthAmerica=nichelandarea_NorthAmerica/grasslandarea_NorthAmerica;
nicheprc_Oceania=nichelandarea_Oceania/grasslandarea_Oceania;
nicheprc_Europe=nichelandarea_Europe/grasslandarea_Europe;
nicheprc_world=nichelandarea_world./grasslandarea_world;


