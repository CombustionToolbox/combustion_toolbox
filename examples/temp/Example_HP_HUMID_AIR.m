% -------------------------------------------------------------------------
% EXAMPLE: HP - Humid air
%
% Compute adiabatic temperature and equilibrium composition at constant
% pressure (e.g., 1.01325 bar) for lean to rich CH4-ideal_air mixtures at
% standard conditions, a set of 10 species considered, and excess of air
% of 15%, and a set of relative humidity (0, 20, 40, 60, 80, 100) [%]
%   
% LS == {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'NO', 'OH', 'O', 'H'}
%   
% See wiki or setListspecies method from ChemicalSystem class for more
% predefined sets of species
%
% Results compared with [1]
%
% [1] Sánchez, Y. A. C., Piedrahita, C. A. R., & Serrano, E. G. F. (2014).
%     Influencia del aire húmedo en la combustión del metano. Scientia et
%     technica, 19(4), 364-370.
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update Feb 21 2024
% -------------------------------------------------------------------------

% User inputs
air_excess = 15;              % [%]
humidity_relative = 0:20:100; % [%]
T = [10:2:60] + 273.15;       % [K]
p = 1 * 1.01325;              % [bar]

% Initialization
self = App();
DB_master = self.DB_master;
DB = self.DB;
LS = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'NO', 'OH', 'O', 'H'};

% Definitions
phi_t = 2; % stoichiometric number of moles of air per mole of fuel
phi = 1 / (1 + air_excess / 100); % equivalence ratio [-]
mm_air = self.DB.Air.mm; % [g/mol]
mm_H2O = self.DB.H2O.mm; % [g/mol]
n_air = 4.76; % total number of moles of stoichiometric ideal dry air (79% N2 and 21% O2 in volume)

% Miscellaneous
config = self.Misc.config;
config.innerposition = [0.2 0.2 0.35 0.5];  % Set figure inner position [normalized]
config.outerposition = [0.2 0.2 0.35 0.5];  % Set figure outer position [normalized]
colors = brewermap(length(humidity_relative), 'Greys');
FLAG_FIRST = true;

for j = length(humidity_relative):-1:1
    % Get specific humidity
    W = humidity_specific(T, p, humidity_relative(j));
    % Get moles H2O
    n_H20 = W * n_air * (mm_air / mm_H2O);
    % Solve problem
    for i = length(T):-1:1
        % Initialization
        self = App('fast', DB_master, DB, LS);
        % Initial conditions
        self = set_prop(self, 'TR', T(i), 'pR', p, 'phi', phi);
        self.PD.S_Fuel     = {'CH4'};
        self.PD.S_Oxidizer = {'N2', 'O2', 'H2O'};
        self.PD.N_Oxidizer = phi_t ./ phi .* [79/21, 1, n_H20(i)];
        % self.PD.ratio_oxidizers_O2 = [79, 21, nH20(i) * 21] / 21;
        % Set additional inputs (depends of the problem selected)
        self = set_prop(self, 'pP', self.PD.pR.value);
        % Solve problem
        self = solve_problem(self, 'HP');
        % Get data
        mix1{i} = self.PS.strR{1};
        mix2{i} = self.PS.strP{1};
    end
    
    if FLAG_FIRST
        index_CO2 = find_ind(self.S.LS, 'CO2');
        index_H2O = find_ind(self.S.LS, 'H2O');
        index_CH4 = find_ind(self.S.LS, 'CH4');
        FLAG_FIRST = false;
    end

    % Extract data
    mR = cell2vector(mix1, 'mi');
    TR = cell2vector(mix1, 'T');
    Yi_R = cell2vector(mix1, 'Yi');
    Yi_R_CH4 = Yi_R(index_CH4, :);

    mP = cell2vector(mix2, 'mi');
    TP = cell2vector(mix2, 'T');
    Yi_P = cell2vector(mix2, 'Yi');
    Xi_P = cell2vector(mix2, 'Xi');
    Yi_P_H2O = Yi_P(index_H2O, :);
    Xi_P_CO2 = Xi_P(index_CO2, :);
    
    m_CH4 = mR .* Yi_R_CH4;
    m_H20 = mP .* Yi_P_H2O;
    
    % Plot results
    if ~exist('ax1', 'var')
        ax1 = set_figure(config);
        ax2 = set_figure(config);
        ax3 = set_figure(config);
    end
    
    % Adiabatic flame temperature vs. reactant's temperature
    ax1 = plot_figure('T_R\ [\rm K]', TR, 'T_P\ [\rm K]', TP, 'ax', ax1, 'color', colors(j, :));
    % Molar fraction of CO2 vs. reactant's temperature
    ax2 = plot_figure('T_R\ [\rm K]', TR, 'X_{\rm CO_2}', Xi_P_CO2, 'ax', ax2, 'color', colors(j, :));
    % Water content after the combustion process (kg H2O / kg fuel) vs. reactant's temperature
    ax3 = plot_figure('T_R\ [\rm K]', TR, '\rm Water\ content\ [kg_{H_2O}/kg_{CH_4}]', m_H20./m_CH4, 'ax', ax3, 'color', colors(j, :));
end