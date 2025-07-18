%% === Step 1 – Extract data ===
monthLabels = string(benchmark_data{:,1});
dataTable    = benchmark_data(:, 2:end);
varNames     = dataTable.Properties.VariableNames;

%% === Step 2 – Select rows & columns ===
selectedRows = 1:12;
selectedCols = [3, 8, 16, 21, 63, 64, 65, 69, 70, 71, 57, 58, 59] - 1;
selectedData = dataTable{selectedRows, selectedCols};
selectedMonthLabels = monthLabels(selectedRows);
selectedVarNames    = varNames(selectedCols);

%% STEP 3 – column labels (safe La TeX)
customColumnLabels = { ...
    'Perfect Foresight', ...
    'MPC-V, $N_{ext}=48$', ...
    'MPC-V, $N_{ext}=72$', ...
    'MPC-V, $N_{ext}=96$', ...
    'MPC-1, $N_{ext}=48$', ...
    'MPC-1, $N_{ext}=72$', ...
    'MPC-1, $N_{ext}=96$', ...
    'MPC-2, $N_{ext}=48$', ...
    'MPC-2, $N_{ext}=72$', ...
    'MPC-2, $N_{ext}=96$', ...
    'MPC-3, $N_{ext}=48$', ...
    'MPC-3, $N_{ext}=72$', ...
    'MPC-3, $N_{ext}=96$'};

%% STEP 4 – create the heat-map (note the new argument)
% --- create the heat-map -------------------------------------------------
h = heatmap( customColumnLabels, ...
             selectedMonthLabels, ...
             selectedData/10, ...
             'Interpreter',      'latex', ...
             'CellLabelFormat', '\\bf %.2f', ...
             'CellLabelColor',   'black', ...   % ⇦ force all labels to black
             'FontColor',        'black');      %   axis ticks / titles

%% === Step 5 – Add titles/labels with LaTeX interpreter ===
% t  = title (h, '$\bf Benchmark\;Data\;with\;Uncertainty$', 'Interpreter','latex');
% xl = xlabel(h, '$\textbf{Type}$',                          'Interpreter','latex');
% yl = ylabel(h, '$\textbf{Month}$',                        'Interpreter','latex');

%% === Step 6 – Tell the underlying axes to interpret tick labels as LaTeX ===
% (The heatmap stores its categorical axes in hidden children.)
ax = struct(h).Axes;
set(ax, 'TickLabelInterpreter','latex');

ax = gca;
ax.FontSize   = 18;
% ax.Box        = 'off';
% ax.XColor     = [0.3 0.3 0.3];
% ax.YColor     = [0.3 0.3 0.3];
%% === Step 7 – Optional appearance tweaks ===
h.ColorbarVisible = 'on';
