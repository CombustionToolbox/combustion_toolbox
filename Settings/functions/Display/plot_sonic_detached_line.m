x_sonic = [0, theta_sonic];
y_sonic = [90, beta_sonic];
x_max = [0, theta_max];
y_max = [90, beta_max];

start1 = 0;
start2 = 35;

degree_2 = 2;

ind_2 = find(x_max > start2, 1, "first");
p_max_2 = polyfit(x_max(ind_2:end), y_max(ind_2:end), degree_2);
p_sonic_2 = polyfit(x_sonic(ind_2:end), y_sonic(ind_2:end), degree_2);
x_max_new_2 = linspace(start2, max(x_max));
y_max_new_2 = polyval(p_max_2, x_max_new_2);
x_sonic_new_2 = linspace(start2, max(x_sonic));
y_sonic_new_2 = polyval(p_sonic_2, x_sonic_new_2);


x_max_new_1 = linspace(0, start2);
y_max_new_1 = interp1(x_max(1:ind_2-1), y_max(1:ind_2-1), x_max_new_1, 'linear');
x_sonic_new_1 = linspace(0, start2);
y_sonic_new_1 = interp1(x_sonic(1:ind_2-1), y_sonic(1:ind_2-1), x_sonic_new_1, 'linear');


x_max_new = [x_max_new_1, x_max_new_2];
y_max_new = [y_max_new_1, y_max_new_2];

x_sonic_new = [x_sonic_new_1, x_sonic_new_2];
y_sonic_new = [y_sonic_new_1, y_sonic_new_2];

y_max_new = spline(x_max_new, y_max_new, 100);
y_sonic_new = spline(x_sonic_new, y_sonic_new, 100);

nfrec = 2;
y_max_new_smooth = smooth(x_max_new(1:nfrec:end), y_max_new(1:nfrec:end));
y_sonic_new= smooth(x_sonic_new, y_sonic_new);

plot(x_max_new(1:nfrec:end), y_max_new_smooth, 'k:', 'LineWidth', self.Misc.config.linewidth);
plot(x_sonic_new(1:nfrec:end), y_sonic_new(1:nfrec:end), 'k:', 'LineWidth', self.Misc.config.linewidth);

plot(x_max, y_max, 'LineWidth', self.Misc.config.linewidth);
plot(x_sonic, y_sonic, 'LineWidth', self.Misc.config.linewidth);