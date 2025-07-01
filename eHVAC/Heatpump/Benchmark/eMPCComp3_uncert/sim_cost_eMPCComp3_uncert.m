clearvars -except meta_data
% addpath('MetaData')
% load('MetaData\meta_data.mat','meta_data')
load('Heatpump\Benchmark\benchmark_data.mat')

Nh_add = 48;

% Loop over months
for year = 2023:2023
    for month = 7:12
        % Get first and last day of the month
        first_day = datetime(year, month, 1);
        last_day = dateshift(first_day, 'end', 'month');
        num_days = day(last_day);

        % Format 'MM-YYYY'
        month_str = datestr(first_day, 'mm-yyyy');

        % Find start index in meta_data
        start_idx = find(meta_data.DateTime == first_day, 1);

        % Run simulation and collect bill
        fprintf('Processing %s...\n', month_str);

        bill = cost_eMPCComp3_uncert(Nh_add, start_idx, month_str, num_days, meta_data);
        bill = mean(bill)

        row_idx = find(benchmark_data.Month == month_str);
        coll_name = ['eMPCComp3_uncert', num2str(Nh_add), 'h_add'];
        benchmark_data.(coll_name)(row_idx) = bill;
        
    end
end


% Save table
save('Heatpump/Benchmark/benchmark_data.mat', 'benchmark_data');

fprintf('✅ Benchmark table saved with %d entries.\n', height(benchmark_data));