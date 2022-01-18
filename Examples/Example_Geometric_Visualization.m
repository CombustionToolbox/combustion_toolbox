% -------------------------------------------------------------------------
% EXAMPLE: Geometric Visualization
%
% Compute equilibrium composition at defined temperature (e.g., 5000 K) and
% pressure (e.g., 1.01325 bar) for a diatomic mixture at standard
% conditions, and a set of 2 species considered.
%   
% species == {'N2', 'N'}
%   
% See wiki or ListSpecies() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                  
% Last update Dec 10 2021
% -------------------------------------------------------------------------

%% INITIALIZE
self = App({'H', 'H2'});
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325);
self.PD.S_Fuel = {'H2', 'H'};
self.PD.N_Fuel = [1, 1];
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'pP', self.PD.pR.value, 'TP', 4000); 
%% SOLVE PROBLEM
self = SolveProblem(self, 'TP');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
%% GEOMETRICAL VISUALIZATION
moles_CT = self.PS.strP{1}.N .* self.PS.strP{1}.Xi;
table(self.S.LS', moles_CT)
clear moleFractions
N = 200;
levels = 5;
x1 = linspace(0, 3, N);
x2 = linspace(0, 2.25, N);
[X1, X2] = meshgrid(x1, x2);
moleFractions(:, :, 1) = X1;
moleFractions(:, :, 2) = X2;
G = compute_G(self, moleFractions);
contour(moleFractions(:, :, 2), moleFractions(:, :, 1), G, 1e4 * [-6, -7.28, -10.433, -12])
hold on;


restriction_fun = sum(self.C.A0.value(:, 1:end)) - 2 * x2;
ind_zero = restriction_fun > 0;
plot(x2(ind_zero), restriction_fun(ind_zero));
plot(moles_CT(2), moles_CT(1), 'ko', 'MarkerFaceColor', 'k');
ylim([0, max(x1)]);
xlim([0, max(x2)])

function G = compute_G(self, moles0)
    [Nx, Ny, Nz] = size(moles0);
    mix2 = self.PS.strP{1};
    isgas = ~self.C.N0.value(:, 2);
    R0TP = self.C.R0 * temperature(mix2); % [J/(mol)]
    g0 = set_g0(self.S.LS, temperature(mix2), self.DB);
    chemPotentials = g0;
    for i = Nx:-1:1
        for j = Ny:-1:1
            moles(:, 1, 1) = moles0(i, j, :);
            chemPotentials_0 = g0(isgas);
            chemPotentials(isgas) = R0TP * log(moles .* isgas / sum(moles)) + log(pressure(mix2));
            G(i, j) = dot(moles, chemPotentials) .* sum(moles);
        end
    end
end