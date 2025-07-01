clearvars -except meta_data
% addpath('MetaData')
% load('MetaData\meta_data.mat','meta_data')
load('Heatpump\Benchmark\OptimalCost\opti_cost_HPeMPC5_04-2023_FULL.mat');


% Generate date vector
date_vec = meta_data.DateTime(start : start+num_days*24-1);
date_vec = repelem(date_vec, H_split);

% Define custom colors
orange = [1 0.5 0];
red    = [1 0 0];

% Clear and create figure
clf
figure(2);

% --- Subplot 1: Temperature States ---
subplot(2,1,1);
hold on;

% Plot temperature-related states
plot(date_vec, X_sim(1,1:end-1), 'b', 'LineWidth', 1.5);              % T_r (Blue)
plot(date_vec, X_sim(2,1:end-1), 'Color', orange, 'LineWidth', 1.5);  % T_f (Orange)
plot(date_vec, X_sim(3,1:end-1), 'Color', red, 'LineWidth', 1.5);     % T_w (Red)
plot(date_vec, Ta(1:Steps), 'g', 'LineWidth', 1);                     % T_a (Green)
plot(date_vec, Tmin_air(1:Steps), 'k--', 'LineWidth', 1);            % T_min (Black dashed)
plot(date_vec, Tmax_air(1:Steps), 'k--', 'LineWidth', 1);            % T_max (Black dashed)

% Legend and labels
legend({'$T_r$', '$T_f$', '$T_w$', '$T_a$', '$T_{\min}$', '$T_{\max}$'}, ...
       'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');
ylabel('Temperature ($^\circ$C)', 'Interpreter', 'latex');
xlim([date_vec(1), date_vec(end)]);
grid on

% Axis formatting
ax = gca;
ax.FontSize   = 16;
ax.LineWidth  = 1.2;
ax.Box        = 'off';
ax.XColor     = [0.3 0.3 0.3];
ax.YColor     = [0.3 0.3 0.3];

% --- Subplot 2: Control Input & Electricity Price ---
subplot(2,1,2);

yyaxis left
stairs(date_vec, U_sim, 'k', 'LineWidth', 1.2);
ylabel('Control Input (W)', 'Interpreter', 'latex');
ax = gca;
ax.YAxis(1).Color = [0 0 0];

yyaxis right
stairs(date_vec, c(1:Steps)*1e-3, 'r', 'LineWidth', 1.5);
ylabel('Electricity Price (€/kWh)', 'Interpreter', 'latex');
ax.YAxis(2).Color = red;

xlim([date_vec(1), date_vec(end)]);
xlabel('Time', 'Interpreter', 'latex');
grid on

% Axis formatting
ax.FontSize   = 16;
ax.LineWidth  = 1.2;
ax.Box        = 'off';
ax.XColor     = [0.3 0.3 0.3];
