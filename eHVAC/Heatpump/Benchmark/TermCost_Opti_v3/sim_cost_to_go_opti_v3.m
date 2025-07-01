clear all
addpath('MetaData')
load('MetaData\meta_data.mat','meta_data')

dt_col = meta_data.DateTime;
tic
% Loop over years and months
for year = 2023:2023
    for month = 6:9

        % Get first and last datetime for this month
        first_day = datetime(year, month, 1);
        last_day = dateshift(first_day, 'end', 'month');

        % Format as 'MM-YYYY' for input to your function
        month_str = datestr(first_day, 'mm-yyyy');

        % Number of days in the month
        num_days = day(last_day);
        
        % Find the index in meta_data where DateTime == first_day
        start_index = find(dt_col == first_day, 1);

        % Skip if not found (e.g., missing data)
        if isempty(start_index)
            fprintf('Start index not found for %s — skipping.\n', month_str);
            continue
        end

        % Call your function
        fprintf('Running for %s (%d days)...\n', month_str, num_days);
        cost_to_go_opti_v3(num_days, month_str, start_index);
    end
end
toc