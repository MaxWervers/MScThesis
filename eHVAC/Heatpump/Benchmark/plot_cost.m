% clearvars -except meta_data
% addpath('MetaData')
% load('MetaData\meta_data.mat','meta_data')
% load('Heatpump\Benchmark\eMPC11Cost\cost_HPeMPC11_48h_01-2023_FULL.mat');
% load('Heatpump\Benchmark\eMPC9Cost_v5\cost_HPeMPC9_v5_48h_add_01-2023_FULL.mat');
% Create date vector and expand for finer resolution
date_vec = meta_data.DateTime(start : start+num_days*24-1);
date_vec = repelem(date_vec, H_split);

% Clear current figure and start figure 2
clf
figure(2);

% Temperature plot
subplot(2,1,1);
hold on;

% Define colors
orange = [1 0.5 0];
red    = [1 0 0];

plot(date_vec, X_sim(1,:), 'b', 'LineWidth', 1.5);              % T_r: Room (Blue)
plot(date_vec, X_sim(2,:), 'Color', orange, 'LineWidth', 1.5);  % T_f: Floor (Orange)
plot(date_vec, X_sim(3,:), 'Color', red, 'LineWidth', 1.5);     % T_w: Water (Red)
plot(date_vec, Ta(1:Steps), 'g', 'LineWidth', 1);               % T_a: Ambient (Green)
plot(date_vec, Tmin_air(1:Steps), 'k--', 'LineWidth', 1);       % T_min
plot(date_vec, Tmax_air(1:Steps), 'k--', 'LineWidth', 1);       % T_max

legend({'$T_r$', '$T_f$', '$T_w$', '$T_a$', '$T_{\min}$', '$T_{\max}$'}, ...
       'Interpreter', 'latex', 'Location', 'best', 'Box', 'off');

ylabel('Temperature ($^\circ$C)', 'Interpreter', 'latex');
xlim([date_vec(1), date_vec(end)]);
grid on

ax = gca;
ax.FontSize   = 16;
ax.LineWidth  = 1.2;
ax.Box        = 'off';
ax.XColor     = [0.3 0.3 0.3];
ax.YColor     = [0.3 0.3 0.3];

% Control and price plot
subplot(2,1,2);

% Left Y-axis: Control Input (Black)
yyaxis left
stairs(date_vec(1:end-1), U_sim, 'k', 'LineWidth', 1.2);
ylabel('Control Input (W)', 'Interpreter', 'latex');
ax = gca;
ax.YAxis(1).Color = [0 0 0];  % Set left Y-axis to black

% Right Y-axis: Electricity Price (Red)
yyaxis right
stairs(date_vec(1:end-1), c(1:Steps-1)*1e-3, 'r', 'LineWidth', 1.5);
ylabel('Electricity Price (€/kWh)', 'Interpreter', 'latex');
ax.YAxis(2).Color = [1 0 0];  % Set right Y-axis to red

% xlabel('Date', 'Interpreter', 'latex');
xlim([date_vec(1), date_vec(end)]);
grid on

% Axis styling
ax.FontSize   = 16;
ax.LineWidth  = 1.2;
ax.Box        = 'off';
ax.XColor     = [0.3 0.3 0.3];


