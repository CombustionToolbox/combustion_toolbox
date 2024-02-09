% -------------------------------------------------------------------------
% Check time dependence with number of chemical species
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Feb 08 2024
% -------------------------------------------------------------------------

% Definitions
nfrec = 10;
NS_min = 7;
phi = 0.5:0.1:2;

% Get list of species
self = App;
LS = find_products(self, {'C', 'H', 'O', 'N'}, 'flag_burcat', true);
LS = unique(['CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'Cbgrb', LS], 'stable');

% Preallocation
NS = NS_min:nfrec:length(LS);
time = zeros(size(NS));
id = 0;

for i = NS_min:nfrec:length(LS)
    % Update
    id = id + 1;
    
    % Initialization
    self = App(LS(1:i));

    % Set initial state
    self = set_prop(self, 'TR', 300, 'pR', 1, 'phi', phi);
    self.PD.S_Fuel     = {'CH4'};
    self.PD.S_Oxidizer = {'N2', 'O2'};
    self.PD.ratio_oxidizers_O2 = [79, 21] / 21;

    % Define additiational constrains (depends of the problem selected)
    self = set_prop(self, 'pP', self.PD.pR.value);

    % Tunning properties
    self.TN.tolN = 1e-14;
    self.TN.tol0 = 1e-4;

    % Solve problem
    t0 = tic;
    self = solve_problem(self, 'HP');
    time(id) = toc(t0);
end

time_per_case = time / length(phi);

% Least squares interpolation
[xData, yData] = prepareCurveData( NS, time_per_case );

% Set up fittype and options.
ft = fittype( 'poly1' );
excludedPoints = xData < 1;
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
ax = set_figure;
h = plot( ax, fitresult, xData, yData, excludedPoints );
legend( h, 'Data', 'Least squares interpolation', 'Location', 'NorthEast', 'Interpreter', 'latex' );

% Label axes
xlabel( 'Number chemical species', 'Interpreter', 'latex' );
ylabel( 'Time per case', 'Interpreter', 'latex' );
grid on