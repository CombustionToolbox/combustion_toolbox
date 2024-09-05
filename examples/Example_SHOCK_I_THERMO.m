% -------------------------------------------------------------------------
% EXAMPLE: SHOCK_I_THERMO
%
% Influence of caloric models on jump conditions in normal shocks
%
% This script examines the effects of different caloric models on the jump 
% conditions (changes in temperature, density, and adiabatic index) 
% encountered in normal shock waves. The models under consideration are:
%
%   1. Thermochemical frozen: assumes a calorically perfect gas, where 
%      specific heat values remain constant.
%   
%   2. Frozen: assumes a thermally perfect gas, where specific heat values 
%      vary with temperature.
%   
%   3. Equilibrium: assumes a calorically imperfect gas, where both 
%      thermal and caloric properties vary, accounting for chemical 
%      reactions in equilibrium.
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                  
% Last update Sep 05 2024
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.shockdetonation.ShockSolver
import combustiontoolbox.utils.display.*

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB);

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'N2', 'O2'}, [79/21, 1]);

% Initialize figure
plotConfig = PlotConfig();
plotConfig.innerposition = [0.05, 0.05, 0.45, 0.55];
plotConfig.outerposition = [0.05, 0.05, 0.45, 0.55];
ax1 = setFigure(plotConfig);
ax2 = setFigure(plotConfig);
ax3 = setFigure(plotConfig);

% Calculations
for i = 1:3

    switch i
        case 1
            FLAG_TCHEM_FROZEN = true;
            FLAG_FROZEN = false;
            linestyle = ':';
        case 2
            FLAG_TCHEM_FROZEN = false;
            FLAG_FROZEN = true;
            linestyle = '--';
        case 3
            FLAG_TCHEM_FROZEN = false;
            FLAG_FROZEN = false;
            linestyle = '-';
    end

    % Define properties
    mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1, 'M1', 1:0.1:10);

    % Initialize solver
    solver = ShockSolver('problemType', 'SHOCK_I', 'FLAG_TCHEM_FROZEN', FLAG_TCHEM_FROZEN, 'FLAG_FROZEN', FLAG_FROZEN);

    % Solve problem
    [mixArray1, mixArray2] = solver.solveArray(mixArray1);
    
    % Plots
    ax1 = plotFigure('M1', [mixArray1.mach], 'T_2/T_1', [mixArray2.T] ./ [mixArray1.T], 'linestyle', linestyle, 'color', [0 0 0], 'ax', ax1);
    ax2 = plotFigure('M1', [mixArray1.mach], '\rho_2/\rho_1', [mixArray2.rho] ./ [mixArray1.rho], 'linestyle', linestyle, 'color', [0 0 0], 'ax', ax2);
    ax3 = plotFigure('M1', [mixArray1.mach], '\gamma_s', [mixArray2.gamma_s], 'linestyle', linestyle, 'color', [0 0 0], 'ax', ax3);

    if i ~= 3
        continue
    end

    ax1 = plotFigure('M1', [5 5], 'T_2/T_1', [ax1.YLim(1), ax1.YLim(2)], 'linestyle', '--', 'color', [0.5 0.5 0.5], 'ax', ax1);
    ax2 = plotFigure('M1', [5 5], '\rho_2/\rho_1', [ax2.YLim(1), ax2.YLim(2)], 'linestyle', '--', 'color', [0.5 0.5 0.5], 'ax', ax2);
    ax3 = plotFigure('M1', [5 5], '\gamma_s', [ax3.YLim(1), ax3.YLim(2)], 'linestyle', '--', 'color', [0.5 0.5 0.5], 'ax', ax3);
end