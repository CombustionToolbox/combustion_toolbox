function check_fitting_thermo(LS, range, property1, strThProp)
% CHECK FITTING OF THERMODYNAMICS PROPERTIES FOR A GIVEN SET OF SPECIES LS
% AND TEMPERATURES.

property1 = lower(property1);

switch property1
    case 'cp'
        property1_curve = 'cPcurve';
    case 'cv'
        property1_curve = 'cVcurve';
    case 'det'
        property1_curve = 'DeTcurve';
        property1       = 'DeT';
    case 'dht'
        property1_curve = 'DhTcurve';
        property1       = 'DhT';
    case 'g0'
        property1_curve = 'g0curve';
    case 'h0'
        property1_curve = 'h0curve';
    case 's0'
        property1_curve = 's0curve';
    otherwise
        error('There is not such property on files');
end

NS = length(LS);
j  = 1; % item legend label
L_legend = cell(1, 2*NS);
for i=NS:-1:1
    % DATA NASA POLYNOMIALS
    y(:, i) = strThProp.(LS{i}).(property1);
    x(:, i) = strThProp.(LS{i}).T;
    % DATA CURVE FITTING
    y_fun{i} = @(T) strThProp.(LS{i}).(property1_curve)(T);
    eval_range(:, i) = feval(y_fun{i}, range);
end

for j=1:NS
    % LEGEND LABELS
    L_legend{j}     = strcat(LS{j},' - NASA');
    L_legend{j + 2} = strcat(LS{j},' - fitting');
end

figure;
hold on;
set(gca, 'XScale', 'log')

plot(x, y, 'linewidth', 1.8)
plot(range, eval_range, ':', 'linewidth', 1.8)
xlabel('$T$', 'fontsize', 16, 'interpreter', 'latex')
ylabel(strcat('$', property1,'$'), 'fontsize', 16, 'interpreter', 'latex')

legend(L_legend, 'fontsize', 14, 'location', 'northeastoutside');
end