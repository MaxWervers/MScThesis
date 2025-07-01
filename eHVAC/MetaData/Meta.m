clear all;
%% Description: Optain datetime and temperature data from a text file
% Define filenames
weather_file1 = 'weer_2011-2020.txt';
weather_file2 = 'weer_2021-2030.txt';

% Read first file
weather_opts1 = detectImportOptions(weather_file1, 'FileType', 'text', 'NumHeaderLines', 32);
weather_opts1.VariableNamesLine = 32;
weather_data1 = readtable(weather_file1, weather_opts1);

% Read second file
weather_opts2 = detectImportOptions(weather_file2, 'FileType', 'text', 'NumHeaderLines', 32);
weather_opts2.VariableNamesLine = 32;
weather_data2 = readtable(weather_file2, weather_opts2);

% Convert YYYYMMDD and HH to MATLAB datetime format
weather_data1.DateTime = datetime(num2str(weather_data1.YYYYMMDD), 'InputFormat', 'yyyyMMdd') + hours(weather_data1.HH);
weather_data2.DateTime = datetime(num2str(weather_data2.YYYYMMDD), 'InputFormat', 'yyyyMMdd') + hours(weather_data2.HH);

% Convert temperature from 0.1 °C to °C
weather_data1.T = weather_data1.T / 10;
weather_data2.T = weather_data2.T / 10;

% Keep only necessary columns
weather_data1 = weather_data1(:, {'DateTime', 'T'});
weather_data2 = weather_data2(:, {'DateTime', 'T'});

% Combine both datasets
weather_data = [weather_data1; weather_data2];

% Sort data by DateTime
weather_data = sortrows(weather_data, 'DateTime');


%% Description: Optain datetime and day-ahead electricity prices data from a CVS file IMPORTANT: PRICE OBTAINED IS THE SAME FOR COMMING HOUR

% Define the list of filenames
elec_fileList = { ...
    'MetaData/elec_prices_2015.csv', ...
    'MetaData/elec_prices_2016.csv', ...
    'MetaData/elec_prices_2017.csv', ...
    'MetaData/elec_prices_2018.csv', ...
    'MetaData/elec_prices_2019.csv', ...
    'MetaData/elec_prices_2020.csv', ...
    'MetaData/elec_prices_2021.csv', ...
    'MetaData/elec_prices_2022.csv', ...
    'MetaData/elec_prices_2023.csv', ...
    'MetaData/elec_prices_2024.csv', ...
    'MetaData/elec_prices_2025.csv', ...
};

% Initialize an empty table to store all data
elec_data = table();

% Loop through each file
for i = 1:length(elec_fileList)
    % Get the current filename
    elec_prices_file = elec_fileList{i};
    
    % Detect import options
    elec_opts = detectImportOptions(elec_prices_file, 'FileType', 'text');
    
    % Read the table
    year_elec_data = readtable(elec_prices_file, elec_opts);
    
    % Identify column names dynamically
    dateColName = elec_opts.VariableNames{1};  
    priceColName = elec_opts.VariableNames{end-2};  
    
    % Split the datetime range at " - "
    splitTimes = split(year_elec_data.(dateColName), ' - ');
    
    % Extract the start time (before "-") and convert to datetime
    year_elec_data.DateTime = datetime(splitTimes(:,1), 'InputFormat', 'dd/MM/yyyy HH:mm:ss');
    
    % Keep only necessary columns
    year_elec_data = year_elec_data(:, {'DateTime', priceColName});
    
    % Append to the main table
    elec_data = [elec_data; year_elec_data]; 
end


%% Merge tables to meta table
% Display variable names in both tables
disp(weather_data.Properties.VariableNames);
disp(elec_data.Properties.VariableNames);

% Merge tables based on Datetime using an outer join
meta_data = outerjoin(weather_data, elec_data, 'Keys', 'DateTime', 'MergeKeys', true);

%% Add temperature bounds
% Define lower and upper bounds based on time and weekday conditions
meta_data.T_LowerBound = 16 * ones(height(meta_data), 1); 
meta_data.T_UpperBound = 24 * ones(height(meta_data), 1); 

% Get hour and weekday info
hours = hour(meta_data.DateTime);
weekdays = weekday(meta_data.DateTime); % 1 = Sunday, 2 = Monday, ..., 7 = Saturday

% Define logical condition for weekdays (Monday to Friday) between 8 AM and 5 PM
work_hours = (weekdays >= 2 & weekdays <= 6) & (hours >= 8 & hours < 17);

% Assign values based on conditions
meta_data.T_LowerBound(work_hours) = 18;
meta_data.T_UpperBound(work_hours) = 22;


%% Add Nh column based on DateTime
hours = hour(meta_data.DateTime);
Nh = zeros(height(meta_data), 1); % Initialize Nh column

% Assign Nh values based on the day-ahead releasse times
Nh(hours == 0) = 24;
Nh(hours == 1) = 23;
Nh(hours == 2) = 22;
Nh(hours == 3) = 21;
Nh(hours == 4) = 20;
Nh(hours == 5) = 19;
Nh(hours == 6) = 18;
Nh(hours == 7) = 17;
Nh(hours == 8) = 16;
Nh(hours == 9) = 15;
Nh(hours == 10) = 14;
Nh(hours == 11) = 13;
Nh(hours == 12) = 12;
Nh(hours == 13) = 35;
Nh(hours == 14) = 34;
Nh(hours == 15) = 33;
Nh(hours == 16) = 32;
Nh(hours == 17) = 31;
Nh(hours == 18) = 30;
Nh(hours == 19) = 29;
Nh(hours == 20) = 28;
Nh(hours == 21) = 27;
Nh(hours == 22) = 26;
Nh(hours == 23) = 25;

% Add Nh column to meta_data
meta_data.Nh = Nh;

%% Split data in train = 0 and validation = 1
% Extract year from DateTime column
years = year(meta_data.DateTime);

% Create new column: 1 for training, 0 for validation
meta_data.TrainVal = ~(years == 2023 | years == 2024); % 1 if not 2023/2024, 0 otherwise

%% Add T_pred_1
% Extract date components
meta_data.Year = year(meta_data.DateTime);
meta_data.Month = month(meta_data.DateTime);
meta_data.Day = day(meta_data.DateTime);
meta_data.Hour = hour(meta_data.DateTime);

% Use only training data (TrainVal == 1)
train_data = meta_data(meta_data.TrainVal == 1, :);

% Compute the average temperature for each (Month, Day, Hour) combination using only training data
avgTemps = groupsummary(train_data, {'Month', 'Day', 'Hour'}, 'mean', 'T');

% Rename the averaged temperature column
avgTemps.Properties.VariableNames{'mean_T'} = 'T_pred_1';

% Join the average temperatures back to the original table
meta_data = outerjoin(meta_data, avgTemps, ...
    'Keys', {'Month', 'Day', 'Hour'}, ...
    'MergeKeys', true);

% Sort by DateTime
meta_data = sortrows(meta_data, 'DateTime');

% Clean up unnecessary columns
meta_data.Month = [];
meta_data.Day = [];
meta_data.Hour = [];
meta_data.Year = [];
meta_data.GroupCount = [];

%% Add T_pred_2
% Extract necessary time components
meta_data.Year = year(meta_data.DateTime);
meta_data.Week = week(meta_data.DateTime);  % Get the week number of the year
meta_data.Hour = hour(meta_data.DateTime);  % Get the hour of the day

% Use only training data (TrainVal == 1)
train_data = meta_data(meta_data.TrainVal == 1, :);

% Compute the average temperature for each (Week, Hour) combination using only training data
weeklyAvgTemps = groupsummary(train_data, ["Week", "Hour"], 'mean', 'T');

% Rename the averaged temperature column
weeklyAvgTemps.Properties.VariableNames{'mean_T'} = 'T_pred_2';

% Merge back into meta_data while keeping only DateTime, Temperature, and WeeklyAvgTemperature
meta_data = outerjoin(meta_data, weeklyAvgTemps, ...
    'Keys', {'Week', 'Hour'}, ...
    'MergeKeys', true);

% Sort by DateTime
meta_data = sortrows(meta_data, 'DateTime');

% Remove temporary columns to maintain original structure
meta_data.Year = [];
meta_data.Week = [];
meta_data.Hour = [];
meta_data.GroupCount = [];

%% Ce_pred_1
meta_data.Ce_pred_1 = cell(height(meta_data), 1); 

% Define parameters
numHours = 100;  % Length of prediction vector
daysBack = 30;   % Lookback period in days
clear hours

% Iterate over each row in the table
for i = 1:height(meta_data)
    currentTime = meta_data.DateTime(i); 
    startHistory = currentTime - daysBack;

    % Extract past data from the last 30 days
    pastData = meta_data(meta_data.DateTime >= startHistory & meta_data.DateTime < currentTime, :);

    % Initialize prediction vector
    predictionVector = zeros(1, numHours);

    % Generate predictions for the next 40 hours
    for j = 0:numHours-1
        futureHour = currentTime + hours(j);

        % Find historical prices for the same hour
        matchingHours = hour(pastData.DateTime) == hour(futureHour);

        % Compute average price
        if any(matchingHours)
            predictionVector(j+1) = mean(pastData.Day_aheadPrice_EUR_MWh_(matchingHours));
        else
            predictionVector(j+1) = NaN; % Handle missing data
        end
    end

    % Store the prediction vector
    meta_data.Ce_pred_1{i} = predictionVector;
end


%% TForecast100hr
% Parameters
forecast_horizon = 200;
phi = 1;
sigma = 0.2;
n_forecasts = 30;  % Number of forecast samples

% Initialize new column
n_rows = height(meta_data);
forecast_cell = cell(n_rows, 1);  % Preallocate cell array

% Loop over valid rows where we can do a 100-hour forecast
for i = 1:(n_rows - forecast_horizon)
    % Extract true temps for the next 100 hours starting at time i
    true_segment = meta_data.T(i:i + forecast_horizon - 1)';

    % Preallocate matrix for forecasts: each row is one forecast
    forecasts = zeros(n_forecasts, forecast_horizon);

    for j = 1:n_forecasts
        % Generate random walk error
        error = zeros(1, forecast_horizon);
        for k = 2:forecast_horizon
            error(k) = phi * error(k - 1) + randn * sigma;
        end

        % Forecast = true + error
        forecasts(j, :) = true_segment + error;
    end

    % Store in cell
    forecast_cell{i} = forecasts;
end

% Fill remaining rows with NaN matrices to maintain same table height
for i = (n_rows - forecast_horizon + 1):n_rows
    forecast_cell{i} = NaN(n_forecasts, forecast_horizon);
end

% Add forecast column to table
meta_data.T_Forecast100hr = forecast_cell;
%%
save("MetaData\meta_data.mat", 'meta_data');
save("MetaData\meta_data.mat", 'meta_data', '-v7.3');
