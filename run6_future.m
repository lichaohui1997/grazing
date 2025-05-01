addpath(genpath('/Users/lichaohui/Desktop/calculation/grazingniche'));
addpath(genpath('/Users/lichaohui/Desktop/calculation/command'));

% imagesc(range_rcp26)
% imagesc(data_hurs_rcp26CHEM)
% imagesc(data_hurs_rcp26ESM)
% imagesc(data_hurs_rcp26MRI)
% imagesc(futureland_rcp26)


%% map area conversion
% mapsize=[720,1440]
% 
% R = 6371; % Earth's radius in km
% latitude_diff = 180 / mapsize(1); % in degrees
% longitude_diff = 360 / mapsize(2); % in degrees, but note that it remains constant around the globe
% 
% % Convert longitude difference to radians
% longitude_diff = deg2rad(longitude_diff);
% 
% % Initialize the areas matrix
% areas = zeros(mapsize(1), 1);
% 
% % Compute area for each latitude band
% for i = 1:mapsize(1)
%     % Calculate the latitude bounds for the current band
%     latitude1 = 90 - (i-1) * latitude_diff; % Start from the top and move downwards
%     latitude2 = 90 - i * latitude_diff;
%     
%     % Convert to radians
%     latitude1 = deg2rad(latitude1);
%     latitude2 = deg2rad(latitude2);
%     
%     % Compute the area and store in the matrix
%     areas(i) = R^2 * longitude_diff * (sin(latitude1) - sin(latitude2));
% end
% 
% % Display the areas matrix
% disp(areas);
% 
% worldarea=repmat(areas,[1,mapsize(2)])
% 
% sum(sum(worldarea))

%% download MRIOC CHEM future rcp data
% this has already been downloaded. does not need to be run twice
% %%
% base_url = 'http://esgf-data01.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM-CHEM';
% 
% % Define the scenarios and variables
% scenarios = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
% variables = {'tas', 'sfcWind', 'hurs', 'pr'};
% pattern = '%s_Amon_MIROC-ESM-CHEM_%s_r1i1p1_200601-210012.nc';
% 
% % Loop over scenarios and variables to download files
% for i = 1:length(scenarios)
%     for j = 1:length(variables)
%         % Construct the file URL and filename
%         filename = sprintf(pattern, variables{j}, scenarios{i});
%         url = strcat(base_url, '/', scenarios{i}, '/mon/atmos/Amon/r1i1p1/v20120710/', variables{j}, '/', filename);
% 
%         % Download the file
%         websave(filename, url);
% 
%         % Display the file content (optional)
%         % ncdisp(filename)
%     end
% end

%% downloading the MRIOC-ESM model

% % this section of code does not need to be run twice
% % Base URL for the data server
% base_url = 'http://esgf-data01.diasjp.net/thredds/fileServer/esg_dataroot/cmip5/output1/MIROC/MIROC-ESM';
% 
% % Define the scenarios and variables
% scenarios = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
% variables = {'tas', 'sfcWind', 'hurs', 'pr'};
% version = 'v20120710'; % Version number
% run = 'r1i1p1'; % Run identifier
% 
% % Pattern for the filename
% pattern = '%s_Amon_MIROC-ESM_%s_%s_200601-210012.nc';
% 
% % Loop over scenarios and variables to download files
% for i = 1:length(scenarios)
%     for j = 1:length(variables)
%         % Construct the filename
%         filename = sprintf(pattern, variables{j}, scenarios{i}, run);
% 
%         % Construct the URL
%         url = strcat(base_url, '/', scenarios{i}, '/mon/atmos/Amon/', run, '/', version, '/', variables{j}, '/', filename);
% 
%         % Download the file
%         fprintf('Downloading %s...\n', filename);
%         websave(filename, url);
% 
%         % Optionally display the file content
%         % ncdisp(filename)
%     end
% end

%% using LUH2 futureland data

ncinfo('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp126-2-1-e_gn_2100-2300.nc')
lon_nc = ncread('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp126-2-1-e_gn_2100-2300.nc', 'lon');
lat_nc = ncread('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp126-2-1-e_gn_2100-2300.nc', 'lat');
[lon_nc, lat_nc] = meshgrid(lon_nc, lat_nc);
% range_rcp26_ = ncread('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp126-2-1-e_gn_2100-2300.nc', 'range');
% range_rcp85_ = ncread('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-MAGPIE-ssp585-2-1-e_gn_2100-2300.nc', 'range');
range_rcp60_ = ncread('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-GCAM-ssp460-2-1-f_gn_2015-2100.nc', 'range');
range_rcp45_ = ncread('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-MESSAGE-ssp245-2-1-f_gn_2015-2100.nc', 'range');
range_rcp26_ = ncread('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-IMAGE-ssp126-2-1-f_gn_2015-2100.nc', 'range');
range_rcp85_ = ncread('multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-MAGPIE-ssp585-2-1-f_gn_2015-2100.nc', 'range');

% range_rcp26=range_rcp26_(:,:,1)' 
% range_rcp85=range_rcp85_(:,:,1)'

range_rcp26=range_rcp26_(:,:,86)'
range_rcp85=range_rcp85_(:,:,86)'
range_rcp45=range_rcp45_(:,:,86)'
range_rcp60=range_rcp60_(:,:,86)'


range_rcp26(range_rcp26>1)=0
range_rcp85(range_rcp85>1)=0
range_rcp45(range_rcp45>1)=0
range_rcp60(range_rcp60>1)=0

R=georefcells([-90,90],[-180,180],size(range_rcp26));
[futureland_rcp26,resizeRrange] = georesize(range_rcp26,R,160/720,"bilinear");
[futureland_rcp85,resizeRrange] = georesize(range_rcp85,R,160/720,"bilinear");
[futureland_rcp45,resizeRrange] = georesize(range_rcp45,R,160/720,"bilinear");
[futureland_rcp60,resizeRrange] = georesize(range_rcp60,R,160/720,"bilinear");

save('./grazingniche/matdata/futureland.mat','futureland_rcp26','futureland_rcp85','futureland_rcp45','futureland_rcp60')

%%
load('futureland.mat') 
%%
% 导入rcp气候数据到变量区，处理数据（units变化）、graph数据
% 为每个情境和变量读取并处理数据
%% MRI data
dir_path = '/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/rcp/MRI';
scenarios = {'rcp26', 'rcp60','rcp45','rcp85'};
variables = {'pr', 'tas', 'hurs', 'sfcWind'};


for s = 1:length(scenarios)
    scenario = scenarios{s};
    for v = 1:length(variables)
        var = variables{v};

        % 文件名
        filename = [var '_Amon_MRI-CGCM3_' scenario '_r1i1p1_200601-210012.nc'];

        % 全路径
        full_path = fullfile(dir_path, filename);

        % 读取数据
        temp_data = ncread(full_path, var, [1 1 1129], [320 160 12]);

        % 根据变量处理数据
        if strcmp(var, 'tas') % 年均温度
            temp_data = rot90(mean(temp_data(:,:,end-11:end), 3)) - 273.15;
        elseif strcmp(var, 'pr') % 降水
            temp_data = rot90(sum(temp_data(:,:,end-11:end), 3)) * 86400 * 30;
        else
            temp_data = rot90(mean(temp_data(:,:,end-11:end), 3));
        end

        % 创建特定的变量名
        variable_name = ['data_' var '_' scenario 'MRI'];
        eval([variable_name ' = temp_data;']);
    end
end
%% MIROC CHEM data 
dir_path = '/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/rcp/MIROCCHEM';
scenarios = {'rcp26', 'rcp60','rcp45','rcp85'};
variables = {'pr', 'tas', 'hurs', 'sfcWind'};


for s = 1:length(scenarios)
    scenario = scenarios{s};
    for v = 1:length(variables)
        var = variables{v};

        % 文件名
        filename = [var '_Amon_MIROC-ESM-CHEM_' scenario '_r1i1p1_200601-210012.nc'];

        % 全路径
        full_path = fullfile(dir_path, filename);

        % 读取数据
        temp_data = ncread(full_path, var, [1 1 1129], [128 64 12]);

        % 根据变量处理数据
        if strcmp(var, 'tas') % 年均温度
            temp_data = rot90(mean(temp_data(:,:,end-11:end), 3)) - 273.15;
        elseif strcmp(var, 'pr') % 降水
            temp_data = rot90(sum(temp_data(:,:,end-11:end), 3)) * 86400 * 30;
        else
            temp_data = rot90(mean(temp_data(:,:,end-11:end), 3));
        end

        % 创建特定的变量名
        variable_name = ['data_' var '_' scenario 'CHEM'];
        eval([variable_name ' = temp_data;']);
    end
end

%% MIROC ESM的数据 
dir_path = '/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/rcp/MIROC';
scenarios = {'rcp26', 'rcp60','rcp45','rcp85'};
variables = {'pr', 'tas', 'hurs', 'sfcWind'};


for s = 1:length(scenarios)
    scenario = scenarios{s};
    for v = 1:length(variables)
        var = variables{v};

        % 文件名
        filename = [var '_Amon_MIROC-ESM_' scenario '_r1i1p1_200601-210012.nc'];

        % 全路径
        full_path = fullfile(dir_path, filename);

        % 读取数据（仅最后12层）
        temp_data = ncread(full_path, var, [1 1 1129], [128 64 12]);

        % 根据变量处理数据
        if strcmp(var, 'tas') % 年均温度
            temp_data = rot90(mean(temp_data(:,:,end-11:end), 3)) - 273.15;
        elseif strcmp(var, 'pr') % 降水
            temp_data = rot90(sum(temp_data(:,:,end-11:end), 3)) * 86400 * 30;
        else
            temp_data = rot90(mean(temp_data(:,:,end-11:end), 3));
        end

        % 创建特定的变量名
        variable_name = ['data_' var '_' scenario 'ESM'];
        eval([variable_name ' = temp_data;']);
    end
end
%% resizing the CHEM and ESM data to be the same as MRI data
% example code:
R = georefcells([-90,90],[-180,180],size(data_hurs_rcp26CHEM));

% [data_hurs_rcp26CHEM, ~] = georesize(data_hurs_rcp26CHEM, R, 160/64, "bilinear");

% Define scenarios, variables, and model types
scenarios = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
variables = {'pr', 'tas', 'hurs', 'sfcWind'};
model_types = {'CHEM', 'ESM'};

% Loop over each scenario, variable, and model type
for i = 1:length(scenarios)
    for j = 1:length(variables)
        for k = 1:length(model_types)
            % Construct the variable name for the data
            data_var_name = sprintf('data_%s_%s%s', variables{j}, scenarios{i}, model_types{k});

            % Check if the data variable exists in the workspace
            if evalin('base', ['exist(''', data_var_name, ''', ''var'')'])
                % Command to resize the data
                resize_cmd = sprintf('[%s, ~] = georesize(%s, R, 160/64, "bilinear");', data_var_name, data_var_name);

                % Evaluate the command in the base workspace
                evalin('base', resize_cmd);
            else
                fprintf('Data variable %s not found in workspace.\n', data_var_name);
            end
        end
    end
end

%% Define the scenarios and variables
scenarios = {'rcp26', 'rcp45', 'rcp60', 'rcp85'};
%scenarios = {'rcp26', 'rcp85'};
variables = {'pr', 'tas', 'hurs', 'sfcWind'};
%variables = {'tas'};

% Loop over each scenario and variable
for i = 1:length(scenarios)
    for j = 1:length(variables)
        % Construct the variable names for the different models
        data_var_MRI = sprintf('data_%s_%sMRI', variables{j}, scenarios{i});
        data_var_CHEM = sprintf('data_%s_%sCHEM', variables{j}, scenarios{i});
        data_var_ESM = sprintf('data_%s_%sESM', variables{j}, scenarios{i});
        data_var_avg = sprintf('data_%s_%s', variables{j}, scenarios{i});

        % Check if the data variables exist in the workspace
        if evalin('base', ['exist(''', data_var_MRI, ''', ''var'')']) && ...
           evalin('base', ['exist(''', data_var_CHEM, ''', ''var'')']) && ...
           evalin('base', ['exist(''', data_var_ESM, ''', ''var'')'])
            % Command to calculate the average
            avg_cmd = sprintf('%s = (%s + %s + %s) ./ 3;', data_var_avg, data_var_MRI, data_var_CHEM, data_var_ESM);

            % Evaluate the command in the base workspace
            evalin('base', avg_cmd);
        else
            fprintf('One or more data variables for %s %s not found in workspace.\n', variables{j}, scenarios{i});
        end
    end
end

%% climate data graph

var_titles = {'Precipitation', 'Mean Annual Temperature', 'Near Surface Wind Speed','Relative Humidity'};

% 输出路径
output_path = '/Users/lichaohui/Desktop/calculation/grazingniche/figures/';
data=data_hurs_rcp26;
% 定义xy轴坐标的刻度
y_ticks = [1, size(data, 1)/2, size(data, 1)];
y_tick_labels = {'90°N', '0°', '90°S'};

x_ticks = [1, size(data, 2)/2, size(data, 2)];
x_tick_labels = {'0°', '180°E', '180°W'};

for s = 1:length(scenarios)
    scenario = scenarios{s};

    % 创建一个新图形窗口
    figure('WindowState','fullscreen');

    % 创建整体标题
    sgtitle(['Scenario: ' scenario], 'FontSize', 16);

    for v = 1:length(variables)
        var = variables{v};

        % 获取特定的变量名
        variable_name = ['data_' var '_' scenario];

        % 从工作空间获取数据
        data = eval(variable_name);

        %for the precipitation data, I want the values to only show for 100-2000 the values over 2000 show the same as 2000
        % If the current variable is Precipitation, clip the values over 2000 to be 2000
if contains(var, 'pr')
    data(data > 2000) = 2000;
end

        % 创建子图
        subplot(2, 2, v);

        % 显示图像
% 显示图像
R = georefcells([-90,90],[0,360],size(data));
ax = worldmap('world');      
setm(ax, 'FFaceColor');    

if contains(var, 'pr')
    custom_colormap = [addcolorplus(338)];
elseif contains(var, 'tas')
    custom_colormap = [addcolorplus(275)];
elseif contains(var, 'hurs')
    custom_colormap = [addcolorplus(288)];
elseif contains(var, 'sfcWind')
    custom_colormap = [addcolorplus(306)];
end
colormap(gca, custom_colormap);  % Set the colormap only for the current subplot

colorbar;
geoshow(flipud(data), R, 'DisplayType', 'texturemap');   
load coastlines;
plotm(coastlat, coastlon, 'Color', 'black'); 


        % 创建子图标题
        title(var_titles{v}, 'FontSize', 14);

        % 设置xy轴刻度
        set(gca, 'YTick', y_ticks, 'YTickLabel', y_tick_labels);
        set(gca, 'XTick', x_ticks, 'XTickLabel', x_tick_labels);

    end

    % 保存图像
    saveas(gcf, fullfile(output_path, ['figure2100_' scenario '.png']));
end
%% future temperature
sgtitle('future temperature', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(316)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
c = colorbar;  
caxis([-50,40])
%c.Ruler.TickLabelFormat='%g%%';
R = georefcells([-90,90],[0,360],size(data_tas_rcp85));
geoshow(flipud(data_tas_rcp85), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
set(gcf, 'Position',  [584,449,684,537])

save("./grazingniche/matdata/future_climate.mat","data_*")
%% 决定futuer niche 并graph niche
load("future_climate.mat")
load("futureland.mat")
% Define the list of scenarios
scenarios = {'rcp26', 'rcp60','rcp45','rcp85'};

% Initialize a structure to store the niche maps
futureniche_struct = struct();

% Loop through each scenario
for s = 1:length(scenarios)
    % Get the scenario
    scenario = scenarios{s};

    % Get the data for the scenario
    data_pr = eval(['data_pr_' scenario]);
    data_tas = eval(['data_tas_' scenario]);
    data_sfcWind = eval(['data_sfcWind_' scenario]);
    data_hurs = eval(['data_hurs_' scenario]);
    data_landuse = eval(['futureland_' scenario]);

% 调整数据
    halfSize = size(data_landuse, 2) / 2;
    data_landuse_new = circshift(data_landuse, [0, halfSize]);
 
        % Define the conditions
    %cond_pr = (data_pr >= 433 & data_pr <= 2368);
    cond_pr = (data_pr >= 50 & data_pr <= 2627);
    cond_tas = (data_tas >= -3 & data_tas <= 29);
    cond_sfcWind = (data_sfcWind >= 1 & data_sfcWind <= 6);
    cond_hurs = (data_hurs >= 39 & data_hurs <= 67);
    cond_landuse = (data_landuse_new >= 0.01);


    % Combine the conditions to create the niche map
    all_conditions = cond_pr & cond_tas & cond_sfcWind & cond_hurs & cond_landuse;


    futureniche = data_landuse_new .* all_conditions;

    % Store the niche map in the structure
    futureniche_struct.(scenario) = futureniche;

    imagesc(futureniche_struct.rcp26)
    imagesc(data_landuse_new)
    imagesc(data_landuse)
    imagesc(all_conditions)

    % Display the niche map
subplot(2, 2, s);
R=georefcells([-90,90],[0,360],size(futureniche));    
sgtitle('Future livestock niche', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1;addcolorplus(332)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
geoshow(flipud(futureniche), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 
    title(['GN distribution under ' scenario]);
    sgtitle('Future grazing niche under scenarios' )
end

set(gcf, 'Position',  [259,455,1140,528])

%saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figuresnew/futureniche_thinnest.svg');

save('./grazingniche/matdata/futureniche.mat','futureniche_struct')

%%

% Define custom colormap
%custom_colormap = [1 1 1; addcolorplus(287)];
custom_colormap = [1 1 1; addcolorplus(332)];

% Set global title for all subplots
sgtitle('Future Livestock Niche', 'FontSize', 16);

% Create subplots
variables = {cond_pr, cond_tas, cond_sfcWind, cond_hurs, cond_landuse};
titles = {'Precipitation', 'Temperature', 'Surface Wind', 'Humidity' 'Grassland'};

for s = 1:5
    subplot(2, 3, s);
    R = georefcells([-90, 90], [0, 360], size(variables{s}));    

    ax = worldmap('world');
    setm(ax, 'FFaceColor', [1 1 1]);
    set(gcf, 'Colormap', custom_colormap);
    geoshow(flipud(variables{s}), R, 'DisplayType', 'texturemap');
    load coastlines
    plotm(coastlat, coastlon, 'Color', 'black'); 
    title(titles{s});
end

% Add a global colorbar (or add individually to each subplot)
colorbar

%%
save('./grazingniche/matdata/futureniche_struct.mat','futureniche_struct')