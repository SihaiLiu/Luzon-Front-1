%% Reading Copernicus data
clc; clear
filepath = 'D:\LZ\';  % Set the file path
file1 = [filepath , 'LZ-sst.nc'];  % Specify the netCDF file name
ncinfo(file1)  % Display information about the netCDF file
lon = ncread(file1,'lon');  % Read longitude data
lat = ncread(file1,'lat');  % Read latitude data
sst1 =  ncread(file1,'analysed_sst') - 273.15;  % Read and convert SST data from Kelvin to Celsius
sst1_time = ncread(file1,'time');  % Read time data
% Convert time data to a date-time vector and create a string array of year and month
sstTimeVec1 = datevec(double(sst1_time)/(3600*24)+datenum(1981,01,01));  
sst1_ym = string([num2str(sstTimeVec1(:,1)), num2str(sstTimeVec1(:,2))]);
sst1_ym_unique = unique(sst1_ym,'stable');  % Get unique year-month combinations
ym_length = numel(sst1_ym_unique);  % Count of unique year-month combinations

% Select a specific area for analysis
x1_pf=117.7250; x2_pf=124.2250; y1_pf=17.7750; y2_pf=21.7250;

% Find indices of the selected box area
box_xrange_pf=find((lon-x1_pf).*(lon-x2_pf)<=0);
box_yrange_pf=find((lat-y1_pf).*(lat-y2_pf)<=0);
lon_sub=lon(box_xrange_pf);
lat_sub=lat(box_yrange_pf);
sst1_sub = sst1(box_xrange_pf, box_yrange_pf, :);

gap = 111*0.05;  % Define the distance between data points

% Define Sobel filters for edge detection (gradient calculation)
sobel_x = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
sobel_y = [1, 2, 1; 0, 0, 0; -1, -2, -1];

% Initialize matrices to store monthly SST data and gradient magnitude
sst1_month = nan([numel(lon), numel(lat), numel(sst1_ym_unique)]); 
front1_month = nan([numel(lon), numel(lat), numel(sst1_ym_unique)]);

% Loop through each month, calculate mean SST and gradient magnitude
for item =1:size(sst1_ym_unique)
    item_find = find(sst1_ym_unique(item) == sst1_ym);
    sst_item = squeeze(mean(sst1(:,:,item_find), 3));
    front1_item = nan(size(sst_item));
    for i = 2:size(lon)-1
        for j = 2:size(lat)-1
            T = sst_item(i-1:i+1, j-1:j+1);  % Extract a 3x3 matrix from SST data
            dTdx = sum(sum(1/4 .* sobel_x .* T));  % Calculate horizontal temperature gradient
            dTdy = sum(sum(1/4 .* T .* sobel_y));  % Calculate vertical temperature gradient
            GM = sqrt(dTdx^2 + dTdy^2) / gap;  % Calculate gradient magnitude
            front1_item(i, j) = GM;  % Assign gradient magnitude to the output matrix
        end 
    end
    front1_month(:,:,item) = front1_item;
    sst1_month(:,:,item)  = sst_item;
end

% Calculate TPI index (Temperature Position Index) by day
% Define the range for TPI calculation, total range is 1.5 degree or about 160km,
% because the resolution of sst is 0.05degree. 
gap = 31; 
mid_gap = (gap-1)/2;
num_lat = numel(lat);
num_lon = numel(lon);
TPI_sub = nan([numel((gap+1)/2:num_lon-mid_gap), numel((gap+1)/2:num_lat-mid_gap), size(sst1_month,3)]);
TPI = nan(size(sst1_month));
bias = (gap+1)/2 - 1; % Offset for the initial index
for i = (gap+1)/2:num_lon-mid_gap
    TPI_sub_mid = [];
    for j = (gap+1)/2:num_lat-mid_gap
        disp([num2str(i), num2str(j)]);
        % Calculate TPI by subtracting the mean SST in the surrounding area
        TPI_sub_mid(j-bias,:) = sst1_month(i, j,:) - nanmean(sst1_month(i-mid_gap:i+mid_gap, j-mid_gap:j+mid_gap,:), [1,2]);
    end
    TPI_sub(i-bias, :,:) = TPI_sub_mid;
end
TPI((gap+1)/2:num_lon-mid_gap, (gap+1)/2:num_lat-mid_gap, :) = TPI_sub;

% Calculate winter climatological mean
winter_find = [3:12:ym_length, 4:12:ym_length, 5:12:ym_length];
clima_sst_mean = nanmean(sst1_month(:,:,winter_find), 3);
clima_SST_mean = nanmean(front1_month(:,:,winter_find), 3);
clima_TPI_mean = nanmean(TPI(:,:,winter_find), 3);
TPI_anomaly = TPI - clima_TPI_mean;

%% Plotting average maps for different months
[LON,LAT] = meshgrid(double(lon),double(lat)); % Create a grid of longitude and latitude
for i = 1:12  
    TPI_mean = mean(TPI(:, :, i:12:488), 3);  % Calculate mean TPI for each month
    this_month = mod(i+9, 12);
    if this_month == 0
        this_month = 12;
    end
    % Plotting TPI
    f1 = figure();
    set(f1, 'visible', 'off');  % Set the figure visibility to off
    m_proj('Mercator','lat',[y1_pf y2_pf],'lon',[x1_pf x2_pf]);  % Set map projection

    m_pcolor(LON,LAT,TPI_mean');  % Plot TPI using pseudocolor plot
    shading interp  % Interpolate colors
    m_coast('patch',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);  % Add coastline
    m_grid('box','on','linest','none','linewidth',1.5,'tickdir','in','backcolor',[1 1 1],'fontname','Time New Roman','fontsize',15)  % Add grid
    colormap(jet)  % Set colormap
    h=colorbar;  % Add colorbar
    set(get(h,'title'),'string','℃');  % Set colorbar title
    title(strcat('TPI-month ', num2str(this_month)))  % Set plot title

    fileName = strcat(filepath, 'TPI\month_average\','LZ-TPI-month ', num2str(this_month), '.tif');  % Specify file name
    saveas(f1, fileName);  % Save the figure
    close(f1);  % Close the figure
end

%% Analysis and plotting for El Nino and La Nina events
% Find indices for El Nino events
elnino_find = find(sst1_ym_unique == "199712" | sst1_ym_unique == "1998 1" | sst1_ym_unique == "1998 2"|...
                    sst1_ym_unique == "200212" | sst1_ym_unique == "2003 1" | sst1_ym_unique == "2003 2"|...
                    sst1_ym_unique == "200912" | sst1_ym_unique == "2010 1" | sst1_ym_unique == "2010 2"|...
                    sst1_ym_unique == "201512" | sst1_ym_unique == "2016 1" | sst1_ym_unique == "2016 2");
% Find indices for La Nina events
lanina_find = find(sst1_ym_unique == "199812" | sst1_ym_unique == "1999 1" | sst1_ym_unique == "1999 2"|...
                    sst1_ym_unique == "200712" | sst1_ym_unique == "2008 1" | sst1_ym_unique == "2008 2"|...
                    sst1_ym_unique == "201012" | sst1_ym_unique == "2011 1" | sst1_ym_unique == "2011 2"|...
                    sst1_ym_unique == "201712" | sst1_ym_unique == "2018 1" | sst1_ym_unique == "2018 2"|...
                    sst1_ym_unique == "202012" | sst1_ym_unique == "2021 1" | sst1_ym_unique == "2021 2");

TPI_elnino = nanmean(TPI_anomaly(:, :, elnino_find), 3);  % Calculate mean TPI for El Nino
TPI_lanina = nanmean(TPI_anomaly(:, :, lanina_find), 3);  % Calculate mean TPI for La Nina

% Perform two-sample paired t-test
[h,p,ci,stats] = ttest(sst_elnino(:), sst_lanina(:), 'Tail', 'left');  % Right tail implies sst_elnino > sst_lanina
[h2,p2,ci2,stats2] = ttest(SST_elnino(:), SST_lanina(:), 'Tail', 'left');
[h3,p3,ci3,stats3] = ttest(TPI_elnino(:), TPI_lanina(:), 'Tail', 'left');

% Plotting TPI for El Nino
f1 = figure();
set(f1, 'visible', 'off');
m_proj('Mercator','lat',[y1_pf y2_pf],'lon',[x1_pf x2_pf]);
m_pcolor(LON,LAT,TPI_elnino');
shading interp
m_coast('patch',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);
m_grid('box','on','linest','none','linewidth',1.5,'tickdir','in','backcolor',[1 1 1],'fontname','Time New Roman','fontsize',15)
colormap(jet)
h=colorbar;
caxis([-0.2 0.2])
set(get(h,'title'),'string','℃');
title('TPI')
fileName = strcat(filepath, 'TPI\enso\','enino-LZ-TPI', '.tif');
saveas(f1, fileName);
close(f1);

% Plotting TPI for La Nina
f1 = figure();
set(f1, 'visible', 'off');
m_proj('Mercator','lat',[y1_pf y2_pf],'lon',[x1_pf x2_pf]);
m_pcolor(LON,LAT,TPI_lanina');
shading interp
m_coast('patch',[.7 .7 .7],'EdgeColor',[.7 .7 .7]);
m_grid('box','on','linest','none','linewidth',1.5,'tickdir','in','backcolor',[1 1 1],'fontname','Time New Roman','fontsize',15)
colormap(jet)
h=colorbar;
caxis([-0.2 0.2])
set(get(h,'title'),'string','℃');
title('TPI')
fileName = strcat(filepath, 'TPI\enso\','lanina-LZ-TPI', '.tif');
saveas(f1, fileName);
close(f1);
