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
%                 
% Last update Aug 05 2024
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.equilibrium.*
import combustiontoolbox.utils.display.*
import combustiontoolbox.utils.findIndex
import combustiontoolbox.utils.extensions.brewermap

% User inputs
airExcess = 15;              % [%]
humidityRelative = 0:20:100; % [%]
T = [10:2:60] + 273.15;      % [K]
p = 1 * 1.01325;             % [bar]

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB, {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'NO', 'OH', 'O', 'H'});

% Initialize mixture
mix = Mixture(system);

% Initialize solver
solver = EquilibriumSolver('problemType', 'HP', 'FLAG_TIME', false);

% Definitions
stoichiometricMoles = 2; % stoichiometric number of moles of air per mole of fuel
equivalenceRatio = 1 / (1 + airExcess / 100); % equivalence ratio [-]
W_air = DB.species.Air.W; % [kg/mol]
W_H2O = DB.species.H2O.W; % [kg/mol]
moles_air = 4.76; % total number of moles of stoichiometric ideal dry air (79% N2 and 21% O2 in volume)

% Miscellaneous
FLAG_FIRST = true;
config = PlotConfig();
config.innerposition = [0.2 0.2 0.35 0.5];  % Set figure inner position [normalized]
config.outerposition = [0.2 0.2 0.35 0.5];  % Set figure outer position [normalized]
colors = brewermap(length(humidityRelative), 'Greys');
ax1 = setFigure(config);
ax2 = setFigure(config);
ax3 = setFigure(config);

for j = length(humidityRelative):-1:1
    % Get specific humidity
    W = humiditySpecific(T, p, humidityRelative(j));

    % Get moles H2O
    moles_H20 = W * moles_air * (W_air / W_H2O);
    
    % Solve problem
    for i = length(T):-1:1
        % Initialize mixture
        mix = Mixture(system);

        % Define chemical state
        set(mix, {'CH4'}, 'fuel', 1);
        set(mix, {'N2', 'O2', 'H2O'}, 'oxidizer', stoichiometricMoles ./ equivalenceRatio .* [79/21, 1, moles_H20(i)]);
        
        % Define properties
        mix1 = setProperties(mix, 'temperature', T(i), 'pressure', p, 'equivalenceRatio', equivalenceRatio);
        mix2 = setProperties(mix, 'temperature', T(i), 'pressure', p, 'equivalenceRatio', equivalenceRatio);

        % Solve problem
        solver.solveArray(mix2);
        
        if FLAG_FIRST
            index_CO2 = findIndex(system.listSpecies, 'CO2');
            index_H2O = findIndex(system.listSpecies, 'H2O');
            index_CH4 = findIndex(system.listSpecies, 'CH4');
            FLAG_FIRST = false;
        end

        % Extract data
        mR(i, j) = mix1.mi;
        TR(i, j) = mix1.T;
        Yi_R_CH4(i, j) = mix1.Yi(index_CH4);
    
        mP(i, j) = mix2.mi;
        TP(i, j) = mix2.T;
        Yi_P_H2O(i, j) = mix2.Yi(index_H2O, :);
        Xi_P_CO2(i, j) = mix2.Xi(index_CO2, :);
        
        m_CH4(i, j) = mR(i, j) .* Yi_R_CH4(i, j);
        m_H20(i, j) = mP(i, j) .* Yi_P_H2O(i, j);
    end
    
    % Plot results
    
    % Adiabatic flame temperature vs. reactant's temperature
    ax1 = plotFigure('T_R\ [\rm K]', TR(:, j), 'T_P\ [\rm K]', TP(:, j), 'ax', ax1, 'color', colors(j, :));
    % Molar fraction of CO2 vs. reactant's temperature
    ax2 = plotFigure('T_R\ [\rm K]', TR(:, j), 'X_{\rm CO_2}', Xi_P_CO2(:, j), 'ax', ax2, 'color', colors(j, :));
    % Water content after the combustion process (kg H2O / kg fuel) vs. reactant's temperature
    ax3 = plotFigure('T_R\ [\rm K]', TR(:, j), '\rm Water\ content\ [kg_{H_2O}/kg_{CH_4}]', m_H20(:, j)./m_CH4(:, j), 'ax', ax3, 'color', colors(j, :));
end