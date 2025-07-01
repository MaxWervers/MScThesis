clearvars -except meta_data

day = 30; 

% The loaded file is assumed to contain a cell array 'results'
file_content = open('cost_to_go_opti_v3_01-2023.mat');
results = file_content.results; 

% Find the cell for which the num_days field equals the desired day.
idx = find(cellfun(@(r) r.num_days, results) == day);

if isempty(idx)
    error('No results found for day = %d.', day);
end

% Use the found result cell (assume only one match)
current_data = results{idx};

bill_values = current_data.bill;    % Array of bill values corresponding to different z_fix values
u_max = current_data.z_fix;       % Corresponding z_fix values

bill_values = bill_values/10;
% u_max = u_max/10;

num_points = length(u_max);
lines = zeros(length(u_max), num_points-1);

for i = 1:num_points-1
    % Compute the slope and intercept of the secant line connecting successive data points.
    slope = (bill_values(i+1) - bill_values(i)) / (u_max(i+1) - u_max(i));
    intercept = bill_values(i) - slope * u_max(i);
    
    % Evaluate the secant line at the points in z_fine.
    lines(:, i) = slope * u_max + intercept;
end

max_approx = max(lines, [], 2);

clf
figure(1);
scatter(u_max, bill_values, 'ro', 'filled', 'DisplayName', 'Terminal Cost');
hold on;

% Plot each secant line (without adding each to the legend)
for i = 1:num_points-1
    plot(u_max, lines(:, i), '--', 'LineWidth', 1, 'HandleVisibility', 'off');
end

% Add ONE representative secant line for the legend
plot(NaN, NaN, '--', 'Color', 'black', 'LineWidth', 1.5, 'DisplayName', 'Secant Lines');

% Plot the final convex approximation (max of secants)
plot(u_max, max_approx, 'b', 'LineWidth', 2, 'DisplayName', 'Max of Secants');



xlabel_handle = xlabel('{\(z(t+N)\)}', 'Interpreter', 'latex');
xlabel_handle.FontSize = 14;

ylabel_handle = ylabel('\bf{Cost (€)}', 'Interpreter', 'tex');
ylabel_handle.FontSize = 14;

% title(sprintf('Convex Approximation Using Max of Linear Secants - %d Days Left', day));

grid on;
ax = gca;
ax.FontSize   = 16;
ax.LineWidth  = 1.2;
ax.Box        = 'off';
ax.XColor     = [0.3 0.3 0.3];
ax.YColor     = [0.3 0.3 0.3];


legend('Interpreter', 'latex','Box','off');
grid on;
hold off;
