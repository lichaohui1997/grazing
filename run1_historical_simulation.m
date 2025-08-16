%% Download climate data
variables = {'pr', 'tas', 'hurs', 'sfcWind'};
%variables = {'huss'};

% start and end year
start_year = 05;
end_year = 80;

% base path
base_path = '/Users/lichaohui/Desktop/calculation/grazingniche/rawdata/climate';

% create weboptions with a larger timeout
options = weboptions('Timeout', 60);

% loop over the variables
for var = variables
    base_url = ['https://ds.nccs.nasa.gov/thredds/fileServer/CMIP5/NASA/GISS/past1000/E2-R_past1000_r1i1p126/' var{:} '_Amon_GISS-E2-R_past1000_r1i1p126_'];
    base_filename = [var{:} '_Amon_GISS-E2-R_past1000_r1i1p126_'];
    for year = start_year:5:(end_year-5)
        % generate the period
        period = sprintf('1%02d101-1%02d012', year, year+5);

        % construct the url and filename
        url = [base_url period '.nc'];
        filename = [base_filename period '.nc'];

        % full path
        full_path = fullfile(base_path, filename);

        % download the file with increased timeout
        %websave(full_path, url, options);
    end

    % Additional loop for the last file with different naming convention
    final_period = sprintf('1%02d101-1%02d012', start_year, start_year+5);
    url = [base_url final_period '.nc'];
    filename = [base_filename final_period '.nc'];
    full_path = fullfile(base_path, filename);
    %websave(full_path, url, options);
end

%% load the variables into directory, graph them, establish geospatial information
% 

for idx = 1:15
    for var = variables
    filename = ['mean_', var{:},'_',num2str(idx), '.mat'];
    load(filename);
    end
end

%% load into directory
% base_path = '/Users/lichaohui/Desktop/calculation/grazingniche/climate';

% Path to save figures
fig_path = '/Users/lichaohui/Desktop/calculation/grazingniche/figures/SI';

% List of variables to process
variables = {'pr', 'tas', 'sfcWind', 'hurs'};

% Loop over all variables
for var = variables
    % Get a list of all .nc files for the current variable
    files = dir(fullfile(base_path, [var{:} '_*.nc']));
    
    % Loop over all files
    for idx = 1:length(files)
        % Full path to the current file
        filename = fullfile(files(idx).folder, files(idx).name);
        
        % Load the data
        data = ncread(filename, var{:}); 
% in this section of code I am reshaping the data into 4 dimensions, so as
% to aggregate the monthly data into yearly data, then calculate the mean
% of first of every 50 years. so the result is pr, hurs, yearly value in year 1500,1550,1600, etc.

% here I am making the above 100 values of relative humidity 100
        if strcmp(var, 'hurs')
            data(data > 100) = 100;
        end


        % Compute the mean for the first time point and convert to Celsius if temperature data
        if strcmp(var, 'tas') % If the variable is 'tas', convert from Kelvin to Celsius
            data = squeeze(mean(reshape(data, 144, 90, 12, []), 3)); % Turn monthly average to yearly average
            mean_data = rot90(data(:,:,1)) - 273.15;
        elseif strcmp(var, 'pr') % If the variable is 'pr', sum up monthly precipitations
            data = squeeze(sum(reshape(data, 144, 90, 12, []), 3)); % Turn monthly aggregate to yearly aggregate
            mean_data = rot90(data(:,:,1));
            mean_data=86400*30*mean_data
        else
            data = squeeze(mean(reshape(data, 144, 90, 12, []), 3)); % Turn monthly average to yearly average
            mean_data = rot90(data(:,:,1));
        end

        % Store the mean_data to a variable
        eval(['mean_' var{:} '_' num2str(idx) ' = mean_data;']);

        % Plot the data
%         figure;
%         imagesc(mean_data);
%         colorbar;
%         title([var{:} ' File ' num2str(idx)]);

        % Save the figure to a PNG file
        saveas(gcf, fullfile(fig_path, ['mean_' var{:} '_' num2str(idx) '.png']), 'png');

        % Save the variable to a .mat file
        save(['mean_' var{:} '_' num2str(idx) '.mat'], ['mean_' var{:} '_' num2str(idx)]);
    end
end

%% load the downloaded historical grazing data into directory and graph the figures
% Path to save figures
fig_path = '/Users/lichaohui/Desktop/calculation/grazingniche/figures/SI';

% Create an array of years
years = 1000:100:1800;

% Loop through the array of years
for i = 1:length(years)
    % Load the data file for the current year
    filename = ['pasture', num2str(years(i)), 'AD.asc'];
    data = importdata(filename);
    
    % Extract the grid values and header information
    gridData = data.data;
    headerInfo = data.textdata;
    
    % Extract header information
    ncols = gridData(1);
    nrows = gridData(2);
    xllcorner = gridData(3);
    yllcorner = gridData(4);
    cellsize = gridData(5);
    
    % Create coordinate vectors
    x = xllcorner + (0:nrows-1) * cellsize;
    y = yllcorner + (0:nrows-1) * cellsize;
    
    % Reshape the grid data into a 2D matrix
    gridMatrix = reshape(gridData(7:size(gridData)), [], nrows)';

    % Display the image
    gridMatrix(gridMatrix < 0) = 0;
    imagesc(gridMatrix);
    colormap([1, 1, 1; jet(255)]); % White color at index 1, followed by the jet colormap
    colorbar;
    
    % Add a title
    title(['Historical Grazing Distribution: Year ', num2str(years(i))]);

    % Save the figure to a PNG file
    saveas(gcf, fullfile(fig_path, ['historical_grazing_distribution_' num2str(years(i)) '.png']), 'png');
    
    % Save the gridMatrix to a variable named gridMatrix and the current year
    assignin('base', ['gridMatrix', num2str(years(i))], gridMatrix);
end

R_hyde = georasterref('LatitudeLimits', [min(y), max(y)], 'LongitudeLimits', [min(x), max(x)], ...
        'RasterSize', [nrows, ncols], 'RasterInterpretation', 'cells');

disp(R_hyde);


%% changing the coordinates of the hyde data from -180-180 to 0-360, 
% Create an array of years
years = 1000:100:1800;

% Loop through the array of years
for i = 1:length(years)
    % Get the matrix for the current year
    gridMatrix = eval(['gridMatrix', num2str(years(i))]);
    
    % Calculate the half size
    halfSize = size(gridMatrix, 2) / 2;
    
    % Shift the matrix
    shiftedMatrix = circshift(gridMatrix, [0, halfSize]);
    
    % Save the shifted matrix to a variable named shiftedMatrix and the current year
    assignin('base', ['shiftedMatrix', num2str(years(i))], shiftedMatrix);
end


imagesc(shiftedMatrix1800);
colormap([1, 1, 1; jet(255)]); % White color at index 1, followed by the jet colormap
colorbar;  
% alright, the graph shows the conversion is successful.



R_hyde = georasterref('LatitudeLimits', [min(y), max(y)], 'LongitudeLimits', [0,360], ...
        'RasterSize', [nrows, ncols], 'RasterInterpretation', 'cells');

disp(R_hyde)

%% Rescale 
numfile=15
resizedpr = cell(numfile,1);
resizedtas = cell(numfile,1);
resizedhurs = cell(numfile,1);
resizedsfcWind = cell(numfile,1);
resizedhyde = cell(numfile,1);

for i = 1:15
    resizedpr{i} = imresize(eval(['mean_pr_', num2str(i)]), [90 180], 'nearest');
    resizedtas{i} = imresize(eval(['mean_tas_', num2str(i)]), [90 180], 'nearest');
    resizedhurs{i} = imresize(eval(['mean_hurs_', num2str(i)]), [90 180], 'nearest');
    resizedsfcWind{i} = imresize(eval(['mean_sfcWind_', num2str(i)]), [90 180], 'nearest');

end

for i = 1:8
    year = 1000 + (i-1) * 100;
    resizedhyde{i} = imresize(eval(['shiftedMatrix', num2str(year)]), [90, 180]);
end

%% historical grazing activity
sgtitle('Historical grazing', 'FontSize', 16);
%custom_colormap = [1 1 1; jet(255)]; 
custom_colormap = [1 1 1; addcolorplus(309)];
ax = worldmap('world')      
setm(ax, 'FFaceColor', [1 1 1]);    
%colormap(ax,custom_colormap);
set(gcf, 'Colormap', custom_colormap);
colorbar
R = georefcells([-90,90],[0,360],size(resizedhyde{8,1}));
geoshow(flipud(resizedhyde{8,1}), R, 'DisplayType', 'texturemap');   
load coastlines
plotm(coastlat, coastlon, 'Color', 'black'); 

saveas(gcf, '/Users/lichaohui/Desktop/calculation/grazingniche/figures/history_grazing.svg');

%% now their coordinates are harmonized, put them in a row vector and draw heat map
resizedpr_vector = reshape(cell2mat(resizedpr(:).'), [],1);
resizedtas_vector = reshape(cell2mat(resizedtas(:).'), [],1);
resizedhurs_vector = reshape(cell2mat(resizedhurs(:).'), [],1);
resizedsfcWind_vector = reshape(cell2mat(resizedsfcWind(:).'), [],1);
resizedhyde_vector = reshape(cell2mat(resizedhyde(:).'), [],1);
resizedpr_vector = reshape(cell2mat(resizedpr(1:2:end).'), [], 1);
resizedtas_vector = reshape(cell2mat(resizedtas(1:2:end).'), [], 1);
resizedhurs_vector = reshape(cell2mat(resizedhurs(1:2:end).'), [], 1);
resizedsfcWind_vector = reshape(cell2mat(resizedsfcWind(1:2:end).'), [], 1);
resizedhyde_vector = reshape(cell2mat(resizedhyde(:).'), [],1);
%% give them standard names so I can juxtapose them with the other climate data
resize_hurs_giss_vector=resizedhurs_vector;
resize_tas_giss_vector=resizedtas_vector;
resize_pr_giss_vector=resizedpr_vector;
resize_sfcWind_giss_vector=resizedsfcWind_vector;
%% save data
resize_hurs_giss_max=resize_hurs_giss_vector
resize_pr_giss_max=resize_pr_giss_vector
resize_tas_giss_max=resize_tas_giss_vector
resize_sfcWind_giss_max=resize_sfcWind_giss_vector
resizedhyde_vector_s=resizedhyde_vector;
resizedhyde_vector_s(resizedhyde_vector_s<=0) = []

save('./grazingniche/matdata/resizedhyde_vector_s.mat','resizedhyde_vector_s')

