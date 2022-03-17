% -------------------------------------------------------------------------
% EXAMPLE: Geometric Visualization
%
% Compute equilibrium composition at defined temperature (e.g., 3500 K) and
% pressure (e.g., 1.01325 bar) for a diatomic mixture at standard
% conditions, and a set of 2 species considered.
%   
% species == {'H2', 'H'}
%   
% See wiki or ListSpecies() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                  
% Last update March 14 2022
% -------------------------------------------------------------------------

%% INITIALIZE
self = App({'H', 'H2'});
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325);
self.PD.S_Fuel = {'H2'};
self.PD.N_Fuel = 1;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
self = set_prop(self, 'pP', self.PD.pR.value, 'TP', 3500); 
%% SOLVE PROBLEM
self = SolveProblem(self, 'TP');
%% DISPLAY RESULTS (PLOTS)
postResults(self);
%% GEOMETRICAL VISUALIZATION
N_CT = self.PS.strP{1}.N;
moles_CT = N_CT * self.PS.strP{1}.Xi;

table(self.S.LS', moles_CT)
clear moles
L = 200;
levels = 50;
moles1 = linspace(0, 2, L);
moles2 = flip(linspace(0, 1, L));
[n1, n2] = meshgrid(moles1, moles2);
moles(:, :, 1) = n1;
moles(:, :, 2) = n2;
mu = compute_chemical_potential(self, moles);

f = figure;
set(f,'units','normalized','innerposition',[0.1 0.1 0.9 0.8],...
    'outerposition',[0.1 0.1 0.9 0.8])
ax = axes(f);
set(ax,'LineWidth', 1.5,'FontSize', 18,'BoxStyle','full')
hold(ax, 'on'); axis(ax, 'tight');

contour3(ax, moles(:, :, 1), moles(:, :, 2), mu, levels)
c = colorbar;
colormap(parula(levels));

total_atoms = sum(self.C.A0.value(:, 1:end) .* [moles1; moles2]);
restriction_fun = (total_atoms - self.C.A0.value(1, 1:end) .* moles1) ./ self.C.A0.value(2, 1:end);

ind_zero = restriction_fun > 0;
plot(moles1(ind_zero), restriction_fun(ind_zero), '--k', 'LineWidth', 1.5);
plot(moles_CT(1), moles_CT(2), 'ko', 'MarkerFaceColor', 'w');
xlim([0, max(moles1)])
ylim([0, max(moles2)]);


xlabel('moles H', 'Interpreter', 'Latex')
ylabel('moles H$_2$', 'Interpreter', 'Latex')

function mu = compute_chemical_potential(self, ni)
    [Nx, Ny, Nz] = size(ni);
    mix2 = self.PS.strP{1};
    isgas = ~self.C.N0.value(:, 2);
    R0TP = self.C.R0 * temperature(mix2); % [J/(mol)]
    g0 = set_g0(self.S.LS, temperature(mix2), self.DB);
    chemPotentials = g0;
    for i = Nx:-1:1
        for j = Ny:-1:1
            moles(:, 1, 1) = ni(i, j, :);
            chemPotentials(isgas) = g0 .* isgas + R0TP * (log(moles .* isgas / sum(moles)) + log(pressure(mix2)));
            mu(i, j) = dot(moles, chemPotentials(isgas));
        end
    end
end