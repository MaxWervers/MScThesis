clearvars -except meta_data

file_content = open('cost_to_go_opti_v3_01-2023.mat');
results = file_content.results; 

clf
figure(1);
hold on;

% To store minima
minima_z = [];
minima_days = [];
minima_bill = [];

for j = 1:length(results)
    r = results{j};
    if isempty(r.u_max)
        continue;
    end

    x = r.z_fix;
    y = r.num_days * ones(size(x));
    z = r.bill;

    % Scatter points with bill-based color
    scatter3(x, y, z, 60, z, 'filled');

    % Plot connecting lines
    plot3(x, y, z, '-', 'Color', [0.5, 0.5, 0.5, 0.5], 'LineWidth', 1.2);

    % Identify minimum bill for this num_days within a tolerance
    min_bill = min(z);
    tolerance = 1e-2;  % Tolerance value, adjust based on your scale
    idx_min = find(z <= min_bill + tolerance); % allows near-minima

    % Among those, pick the one with max z_fix
    [~, idx_max_z] = max(x(idx_min));
    best_idx = idx_min(idx_max_z);

    % Store the values
    minima_z(end+1) = x(best_idx);
    minima_days(end+1) = r.num_days;
    minima_bill(end+1) = z(best_idx);

end

% Highlight minimum points
scatter3(minima_z, minima_days, minima_bill, 90, 'black', 'filled', 'MarkerEdgeColor', 'k');

ax = gca;
ax.FontSize = 16;  % increases tick label size

yt = ax.YTick;  % current Y tick positions
ytl = arrayfun(@num2str, fliplr(yt), 'UniformOutput', false);  % reverse the labels
ax.YTickLabel = ytl;  % apply the reversed labels

ylim([0 31]);
xlabel('$z_{\mathrm{min}}$','Interpreter','latex', 'FontSize',28);
ylabel('$\tau$','Interpreter','latex', 'FontSize',28);
zlabel('Cost (€)','Interpreter','latex', 'FontSize',22);
% title('Terminal cost function for each day - January 2023','FontSize', 14, 'FontWeight', 'bold');
grid on;
view(45-70, 10);
colormap turbo;
colorbar;