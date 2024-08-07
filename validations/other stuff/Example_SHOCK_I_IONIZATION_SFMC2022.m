% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_I_IONIZATION SFMC 2022
%
% Compute pre-shock and post-shock state for a planar incident shock wave
% at standard conditions, a set of 39 species considered and a set of
% initial shock front velocities (u1) contained in (360, 20000) [m/s]
%    
% Air_ions == {'O2','N2','O','O3','N','NO','NO2','NO3','N2O','N2O3',...
%              'N2O4','N3','eminus','Nminus','Nplus','NOplus','NO2minus',...
%              'NO3minus','N2plus','N2minus','N2Oplus','Oplus','Ominus',...
%              'O2plus', 'O2minus,'CO2','CO','COplus','C','Cplus',...
%              'Cminus','CN','CNplus','CNminus','CNN','NCO','NCN','Ar',...
%              'Arplus'}
%   
% See wiki or ListSpecies() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Oct 22 2021
% -------------------------------------------------------------------------

%% INITIALIZE
self = App('Air_ions');
%% INITIAL CONDITIONS
self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325);
self.PD.S_Oxidizer = {'N2', 'O2', 'Ar', 'CO2'};
self.PD.N_Oxidizer = [78.084, 20.9476, 0.9365, 0.0319] ./ 20.9476;
%% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
u1 = logspace(2, 5, 500); u1 = u1(u1<20000); u1 = u1(u1>=360);
self = set_prop(self, 'u1', u1, 'phi', self.PD.phi.value(1) * ones(1, length(u1)));
%% SOLVE PROBLEM
self = solve_problem(self, 'SHOCK_I');
%% DISPLAY RESULTS (PLOTS)
self.Misc.display_species = self.S.LS;
post_results(self);
%%

Xi = cell2vector(self.PS.strP, 'Xi');
M1 = cell2vector(self.PS.strR, 'u') ./ cell2vector(self.PS.strR, 'sound');
R = cell2vector(self.PS.strP, 'rho') ./ cell2vector(self.PS.strR, 'rho');
P = cell2vector(self.PS.strP, 'p') ./ cell2vector(self.PS.strR, 'p');

ind_pass = sum(Xi > 1e-2, 2) > 0;
Xi = Xi(ind_pass, :);
[Nx, Ny] = size(Xi);
dXidM1 = zeros(Nx, Ny-1);


for i = 1:Nx
    dXidM1(i, 1:Ny-1) = compute_first_derivative(Xi(i, :), M1);
end

% ind_pass = sum(dXidM1 > 0.001, 2) > 0;

% xx = linspace(M1_new(1), M1_new(end), 250);
% dXO2dM1 = interp1(M1_new, dXO2dM1, xx);
% dXN2dM1 = interp1(M1_new, dXN2dM1, xx);
% M1_new = xx;

colorbw = brewermap(Nx, self.Misc.config.colorpalette);
ax = set_figure();

for i=1:Nx
    dl = plot_figure('M1', M1(1:end-1), 'dXidM1', dXidM1(i, :), 'color', colorbw(i,:), 'XScale', 'log', 'ax', ax);
end
legend(self.Misc.display_species(ind_pass))