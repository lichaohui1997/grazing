addpath(genpath('/Users/lichaohui/Desktop/calculation/grazingniche'));
addpath(genpath('/Users/lichaohui/Desktop/calculation/command'));

%% Paths
base_dir   = '/Users/lichaohui/Desktop/calculation/grazingniche';
dir_raw    = fullfile(base_dir,'rawdata','climate');
dir_out    = fullfile(base_dir,'matdata');
if ~exist(dir_raw,'dir'); mkdir(dir_raw); end
if ~exist(dir_out,'dir'); mkdir(dir_out); end

% Variables & years
variables = {'pr','tas','sfcWind','hurs'};   % order matters for scaling/offset
years     = 2010:2018;
months    = 1:12;

% Scaling & offsets aligned with 'variables' order
%   pr: divide by 100
%   tas: divide by 10 then add -273.15 (K -> °C)
%   sfcWind: divide by 1000
%   hurs: divide by 100
scaling_factors = [100, 10, 1000, 100];
offsets         = [  0, -273.15,   0,   0];

% Target grid size (1740 + 60-pad = 1800 rows total)
target_rows = 1740; target_cols = 3600;
north_pad   = 60;   % to prepend zeros, making final 1800 x 3600

% Download & process
for vi = 1:numel(variables)
    var   = variables{vi};
    scale = scaling_factors(vi);
    offs  = offsets(vi);

    % accumulate ANNUAL aggregates into a 3D stack: [rows, cols, years]
    annual_stack = [];

    % --------- DOWNLOAD (skip if file exists) ---------
    for yy = years
        for mm = months
            url = sprintf('https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/monthly/%s/CHELSA_%s_%02d_%d_V.2.1.tif', ...
                          var, var, mm, yy);
            fname = sprintf('CHELSA_%s_%02d_%d_V.2.1.tif', var, mm, yy);
            fpath = fullfile(dir_raw, fname);
            if ~isfile(fpath)
                try
                    websave(fpath, url);
                catch ME
                    warning('Failed to download %s: %s', url, ME.message);
                end
            end
        end
    end

    % --------- BUILD ANNUAL AGGREGATES (per year) ---------
    for yy = years
        % Preallocate monthly stack for this year
        monthly_stack = zeros(target_rows, target_cols, numel(months), 'double');

        for mi = 1:numel(months)
            mm = months(mi);
            fname = sprintf('CHELSA_%s_%02d_%d_V.2.1.tif', var, mm, yy);
            fpath = fullfile(dir_raw, fname);
            if ~isfile(fpath)
                error('Missing file: %s', fpath);
            end

            % Load, resize, scale, offset
            raw = double(imread(fpath));
            raw = imresize(raw, [target_rows, target_cols]);  % 1740 x 3600
            monthly_stack(:,:,mi) = raw./scale + offs;
        end

        % Annual aggregate: sum for precipitation; mean for others
        if strcmp(var,'pr')
            annual = sum(monthly_stack, 3, 'omitnan');
        else
            annual = mean(monthly_stack, 3, 'omitnan');
        end

        % Append northern padding to reach 1800 x 3600
        annual = [zeros(north_pad, target_cols); annual];

        % Accumulate across years
        if isempty(annual_stack)
            annual_stack = zeros(size(annual,1), size(annual,2), numel(years));
        end
        annual_stack(:,:,yy - years(1) + 1) = annual;
    end

    % --------- 9-year MEAN (2010–2018) ---------
    nine_year_mean = mean(annual_stack, 3, 'omitnan');  % 1800 x 3600

    % Save
    aggregate_data = nine_year_mean; %#ok<NASGU>
    save(fullfile(dir_out, sprintf('aggregate_%s.mat', var)), 'aggregate_data','-v7.3');

    fprintf('Saved 2010–2018 mean for %s to %s\n', var, fullfile(dir_out, sprintf('aggregate_%s.mat', var)));
end

%% Load the data
% if you are running the code not for the first time, you only need to
% start with this section
load('./grazingniche/matdata/aggregate_hurs.mat')
aggregate_hurs=aggregate_data;
load('./grazingniche/matdata/aggregate_pr.mat')
aggregate_pr=aggregate_data;
load('./grazingniche/matdata/aggregate_sfcWind.mat')
aggregate_sfcWind=aggregate_data;
aggregate_sfcWind(aggregate_sfcWind==66)=0
load('./grazingniche/matdata/aggregate_tas.mat')
aggregate_tas=aggregate_data;

%% grassland data 
% this section does not need to be run twice
% if you are running the code for a second time, please run the below section code
% and not this one
% please see the documentation of the data here: https://www.earthenv.org/landcover
% shrubs=readgeoraster('consensus_full_class_5.tif')
% (save('shrubs.mat', 'shrubs'))
% herbaceous=readgeoraster('consensus_full_class_6.tif')
% (save('herbaceous.mat', 'herbaceous'))
% cultivated=readgeoraster('consensus_full_class_7.tif')
% (save('cultivated.mat', 'cultivated'))
% builtup=readgeoraster('consensus_full_class_9.tif')
% (save('builtup.mat', 'builtup'))
% barren=readgeoraster('consensus_full_class_11.tif')
% (save('barren.mat', 'barren'))
%
% load('shrubs.mat')
% load('herbaceous.mat')
% load('cultivated.mat')
% load('builtup.mat')
% load('barren.mat')
% Here I am processing the new land use data. 
% This dataset only covers 56S-90N. I am adding the south grid cells
% herbaceous=double(herbaceous);
% shrubs=double(shrubs);
% cultivated=double(cultivated);
% builtup=double(builtup);
% barren=double(barren);
% 
% gap=zeros(4080,43200)

% herbaceous_full=[herbaceous;gap];
% save('./grazingniche/matdata/herbaceous_full.mat', 'herbaceous_full', '-v7.3');
% shrubs_full=[shrubs;gap];
% save('./grazingniche/matdata/shrubs_full.mat', 'shrubs_full', '-v7.3');
% cultivated_full=[cultivated;gap];
% save('./grazingniche/matdata/cultivated_full.mat', 'cultivated_full', '-v7.3');
% builtup_full=[builtup;gap];
% save('./grazingniche/matdata/builtup_full.mat', 'builtup_full', '-v7.3');
% barren_full=[barren;gap];
% save('./grazingniche/matdata/barren_full.mat', 'barren_full', '-v7.3');
%% load the new land use data
% Or else you can just download the ('grassland_env.mat') data from below
% section and ignore this section
% or load the data from here and neglect the above section
% this section takes around 30s to run
load('herbaceous_full.mat')
load('shrubs_full.mat')
% in case you need the other land use types
% load('cultivated_full.mat')
%load('builtup_full.mat')
load('barren_full.mat')
%% resize the loaded grassland maps (herbaceous+shrubs)
R_env = georefcells([-90,90],[-180,180],size(herbaceous_full));
[herbaceous_resized,R_resized] = georesize(herbaceous_full,R_env,1/12,"bilinear");
[shrubs_resized,R_resized] = georesize(shrubs_full,R_env,1/12,"bilinear");
[barren_resized,R_resized]=georesize(barren_full,R_env,1/12,"bilinear")
%[builtup_resized,R_resized]=georesize(builtup_full,R_env,1/12,"bilinear")
%% 
imagesc(barren_resized)
colorbar
%% obtaining a grassland map called grassland_env
grassland_env=shrubs_resized+herbaceous_resized;
save('./grazingniche/matdata/grassland_env.mat','grassland_env');

%% Climate Niche
% load('grassland_env.mat')
imagesc(grassland_env)
colorbar
%% grassland mapping
% to run this you need to add the command file into directory
figure1 = figure('WindowState','fullscreen');
sgtitle('Grassland map', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(332)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[-180,180],size(grassland_env));
geoshow(flipud(grassland_env), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

% set(gcf,'renderer','painters');
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/grassland.svg');

% %% examining the niche
% 
% figure;
% t = tiledlayout(2, 2);
% 
% % Plot mean annual temperature
% nexttile;
% imagesc(aggregate_tas);
% colorbar;
% caxis([5 29]);
% title('Mean Annual Temperature');
% 
% % Plot precipitation
% nexttile;
% imagesc(aggregate_pr);
% caxis([324 2805]);
% colorbar;
% title('Precipitation');
% 
% % Plot near surface wind speed
% nexttile;
% imagesc(aggregate_sfcWind);
% colorbar;
% caxis([1.5 6.0]);
% title('Near Surface Wind Speed');
% 
% % Plot relative humidity
% nexttile;
% imagesc(aggregate_hurs);
% colorbar;
% caxis([47 87]);
% title('Relative Humidity');
% 
% % Adjust the spacing and title of the tiled layout
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
% title(t, 'Climate Data');
% 
% % Save the figure as a PNG
% saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/niche_range.png');
%%
% To change the color scheme so that values within a 
% specified range are shown in green and values outside 
% this range are shown in white, you can create a custom 
% colormap and use logical indexing.
figure;
t = tiledlayout(2, 2);

% Custom colormap
colormap_matrix = [1 1 1; 0 1 0]; % White and green

% Plot mean annual temperature
nexttile;
imagesc(aggregate_tas);
colormap(gca, colormap_matrix);
mask = aggregate_tas >= -3 & aggregate_tas <= 29;
caxis([0 1]);
image_alpha = ones(size(mask));
image_alpha(~mask) = 0;
alpha(image_alpha);
colorbar;
title('Mean Annual Temperature');

% Plot precipitation
nexttile;
imagesc(aggregate_pr);
colormap(gca, colormap_matrix);
mask = aggregate_pr >= 50 & aggregate_pr <= 2627;
caxis([0 1]);
image_alpha = ones(size(mask));
image_alpha(~mask) = 0;
alpha(image_alpha);
colorbar;
title('Precipitation');

% Plot near surface wind speed
nexttile;
imagesc(aggregate_sfcWind);
colormap(gca, colormap_matrix);
mask = aggregate_sfcWind >= 1 & aggregate_sfcWind <= 6;
caxis([0 1]);
image_alpha = ones(size(mask));
image_alpha(~mask) = 0;
alpha(image_alpha);
colorbar;
title('Near Surface Wind Speed');

% Plot relative humidity
nexttile;
imagesc(aggregate_hurs);
colormap(gca, colormap_matrix);
mask = aggregate_hurs >= 39 & aggregate_hurs <= 67;
caxis([0 1]);
image_alpha = ones(size(mask));
image_alpha(~mask) = 0;
alpha(image_alpha);
colorbar;
title('Relative Humidity');

% Adjust the spacing and title of the tiled layout
t.TileSpacing = 'compact';
t.Padding = 'compact';
title(t, 'Climate Data');

%% couple the niche with the land use data
% in order to run this section you need to add the command file into
% directory in order for addcolorplus to work

cond1_pr = (aggregate_pr >= 50 & aggregate_pr <= 2627);
cond1_tas = (aggregate_tas >= -3 & aggregate_tas <= 29);
cond1_sfcWind = (aggregate_sfcWind >= 1 & aggregate_sfcWind <= 6);
cond_hurs = (aggregate_hurs >= 39 & aggregate_hurs <= 67);
cond1_landuse = (grassland_env>0);

% Comerbine the conditions to create the niche map
niche1 = cond1_pr & cond1_tas & cond1_sfcWind & cond_hurs & cond1_landuse;

% Create the coupled landuse map
landuse_coupled1 = double(niche1) .* grassland_env;

% producing the graph of grazing niche
R = georefcells([-90,84],[-180,180],size(landuse_coupled1));

figure1 = figure('WindowState','fullscreen');
sgtitle('Niche map', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(332)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
geoshow(flipud(landuse_coupled1), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/niche_general.svg');

save('./grazingniche/matdata/niche.mat',"landuse_coupled1")
