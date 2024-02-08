% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_POLAR_LIMIT_RR_ALTITUDE
%
% Compute limit regular reflections at different altitudes considering a
% thermochemical frozen gas, a chemically frozen gas, and dissociation,
% ionization, vibrational excitation and electronic excitation. The
% calculations are carried out for a set of initial shock front velocities
% u1/a1 = [1.75:0.1:20]
%    
% Air_ions == {'eminus', 'Ar', 'Arplus', 'C', 'Cplus', 'Cminus', ...
%              'CN', 'CNplus', 'CNminus', 'CNN', 'CO', 'COplus', ...
%              'CO2', 'CO2plus', 'C2', 'C2plus', 'C2minus', 'CCN', ...
%              'CNC', 'OCCN', 'C2N2', 'C2O', 'C3', 'C3O2', 'N', ...
%              'Nplus', 'Nminus', 'NCO', 'NO', 'NOplus', 'NO2', ...
%              'NO2minus', 'NO3', 'NO3minus', 'N2', 'N2plus', ...
%              'N2minus', 'NCN', 'N2O', 'N2Oplus', 'N2O3', 'N2O4', ...
%              'N2O5', 'N3', 'O', 'Oplus', 'Ominus', 'O2', 'O2plus', ...
%              'O2minus', 'O3'}
%   
% See wiki or list_species() for more predefined sets of species
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Dec 09 2022
% -------------------------------------------------------------------------

% Definitions
S_Oxidizer = {'N2', 'O2', 'Ar', 'CO2'};
N_Oxidizer = [78.084, 20.9476, 0.9365, 0.0319] ./ 20.9476;
Mach_number = 1.75:0.1:20.5;
z = [0, 15000, 30000];

% Default
FLAG_FROZEN = false;
FLAG_TCHEM_FROZEN = false;

% Calculations
for i = 3:-1:1

    switch i
        case 1
            LS = 'Air_ions';
        case 2
            LS = {'N2', 'O2', 'Ar', 'CO2'};
            FLAG_FROZEN = true;
        case 3
            LS = {'N2', 'O2', 'Ar', 'CO2'};
            FLAG_TCHEM_FROZEN = true;
    end

    for j = length(z):-1:1
        % Initialize
        self = App(LS);
        % Initial conditions
        self.TN.FLAG_FROZEN = FLAG_FROZEN;
        self.TN.FLAG_TCHEM_FROZEN = FLAG_TCHEM_FROZEN;
        [~, ~, TR, pR] = atmos(z(j));
        pR = convert_Pa_to_bar(pR);
        self = set_prop(self, 'TR', TR, 'pR', pR);
        sound_velocity = compute_sound(TR, pR, S_Oxidizer, N_Oxidizer, 'self', self);
        self.PD.S_Oxidizer = S_Oxidizer;
        self.PD.N_Oxidizer = N_Oxidizer;
        % Additional inputs (depends of the problem selected)
        self = set_prop(self, 'u1', sound_velocity * Mach_number);
        % Solve problem
        self = solve_problem(self, 'SHOCK_POLAR_LIMITRR');
        % Get results
        theta_limitRR(:, i, j) = cell2vector(self.PS.str2_1, 'theta');
        beta_limitRR(:, i, j) = cell2vector(self.PS.str2_1, 'beta');
    end

end

% Plots

% Plot settings
% colors = [209 237 179; 172 222 188; 125 203 195; 66 164 200; 20 111 172] / 255;
colors = [175 221 233; 95 188 211; 0 102 128] / 255;

% Plot wave angle [deg] (limit regular reflection) against pre-shock Mach number
ax = set_figure();
for j = 1:length(z)
    ax = plot_figure('Pre-shock Mach number', Mach_number, 'Wave angle limit RR [deg]', beta_limitRR(:, 1, j), 'linestyle', 'k-', 'ax', ax, 'color', colors(j, :));
    ax = plot_figure('Pre-shock Mach number', Mach_number, 'Wave angle limit RR [deg]', beta_limitRR(:, 2, j), 'linestyle', 'k--', 'ax', ax, 'color', colors(j, :));
    ax = plot_figure('Pre-shock Mach number', Mach_number, 'Wave angle limit RR [deg]', beta_limitRR(:, 3, j), 'linestyle', 'k:', 'ax', ax, 'color', colors(j, :));
end

% Plot deflection angle [deg] (limit regular reflection) against pre-shock Mach number
ax = set_figure();
for j = 1:length(z)
    ax = plot_figure('Pre-shock Mach number', Mach_number, 'Deflection angle limit RR [deg]', theta_limitRR(:, 1, j), 'linestyle', 'k-', 'ax', ax, 'color', colors(j, :));
    ax = plot_figure('Pre-shock Mach number', Mach_number, 'Deflection angle limit RR [deg]', theta_limitRR(:, 2, j), 'linestyle', 'k--', 'ax', ax, 'color', colors(j, :));
    ax = plot_figure('Pre-shock Mach number', Mach_number, 'Deflection angle limit RR [deg]', theta_limitRR(:, 3, j), 'linestyle', 'k:', 'ax', ax, 'color', colors(j, :));
end