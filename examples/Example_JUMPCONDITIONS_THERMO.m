% -------------------------------------------------------------------------
% EXAMPLE: JUMPCONDITIONS_THERMO
%
% Influence of caloric models on jump conditions in normal shocks
%
% This script examines the effects of different caloric models on the jump 
% conditions (changes in temperature, density, adiabatic index, and the 
% dimensionless Rankine-Hugoniot slopes) encountered in normal shock waves.
% The models under consideration are:
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
% Last update April 19 2025
% -------------------------------------------------------------------------

% Import packages
import combustiontoolbox.databases.NasaDatabase
import combustiontoolbox.core.*
import combustiontoolbox.shockdetonation.*
import combustiontoolbox.utils.display.*

% Definitions
FLAG_TCHEM_FROZEN = false;
FLAG_FROZEN = false;
FLAG_PAPER = true;

% Get Nasa database
DB = NasaDatabase();

% Define chemical system
system = ChemicalSystem(DB, 'air_ions');

% Initialize mixture
mix = Mixture(system);

% Define chemical state
set(mix, {'N2', 'O2'}, [79, 21] / 21);

% Initialize figure
plotConfig = PlotConfig();
plotConfig.innerposition = [0.05, 0.05, 0.45, 0.55];
plotConfig.outerposition = [0.05, 0.05, 0.45, 0.55];
ax1 = setFigure(plotConfig); ax2 = setFigure(plotConfig);
ax3 = setFigure(plotConfig); ax4 = setFigure(plotConfig);
ax5 = setFigure(plotConfig); ax6 = setFigure(plotConfig);

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
    mixArray1 = setProperties(mix, 'temperature', 300, 'pressure', 1, 'M1', 1:0.05:20);

    % Initialize solver
    solver = JumpConditionsSolver('FLAG_TCHEM_FROZEN', FLAG_TCHEM_FROZEN, 'FLAG_FROZEN', FLAG_FROZEN, 'FLAG_PAPER', FLAG_PAPER, 'FLAG_RESULTS', false);

    % Solve problem
    jumpData = solver.solve(mixArray1);
    
    % Plots
    ax1 = plotFigure('M1', jumpData.M1, 'Tratio', jumpData.Tratio, 'linestyle', linestyle, 'color', [0 0 0], 'ax', ax1);
    ax2 = plotFigure('M1', jumpData.M1, 'Rratio', jumpData.Rratio, 'linestyle', linestyle, 'color', [0 0 0], 'ax', ax2);
    ax3 = plotFigure('M1', jumpData.M1, '\gamma_s', jumpData.gamma2, 'linestyle', linestyle, 'color', [0 0 0], 'ax', ax3);
    ax4 = plotFigure('M1', jumpData.M1, '-Gammas1', -jumpData.Gammas1, 'linestyle', linestyle, 'color', [0 0 0], 'ax', ax4);
    ax5 = plotFigure('M1', jumpData.M1, 'Gammas2', jumpData.Gammas2, 'linestyle', linestyle, 'color', [0 0 0], 'ax', ax5);
    ax6 = plotFigure('M1', jumpData.M1, 'Gammas3', jumpData.Gammas3, 'linestyle', linestyle, 'color', [0 0 0], 'ax', ax6);

    if i ~= 3
        continue
    end

    ax1 = plotFigure('M1', [5 5], 'Tratio', [ax1.YLim(1), ax1.YLim(2)], 'linestyle', '--', 'color', [0.5 0.5 0.5], 'ax', ax1);
    ax2 = plotFigure('M1', [5 5], 'Rratio', [ax2.YLim(1), ax2.YLim(2)], 'linestyle', '--', 'color', [0.5 0.5 0.5], 'ax', ax2);
    ax3 = plotFigure('M1', [5 5], '\gamma_s', [ax3.YLim(1), ax3.YLim(2)], 'linestyle', '--', 'color', [0.5 0.5 0.5], 'ax', ax3);
    ax4 = plotFigure('M1', [5 5], '-Gammas1', [ax4.YLim(1), ax4.YLim(2)], 'linestyle', '--', 'color', [0.5 0.5 0.5], 'ax', ax4);
    ax5 = plotFigure('M1', [5 5], 'Gammas2', [ax5.YLim(1), ax5.YLim(2)], 'linestyle', '--', 'color', [0.5 0.5 0.5], 'ax', ax5);
    ax6 = plotFigure('M1', [5 5], 'Gammas3', [ax6.YLim(1), ax6.YLim(2)], 'linestyle', '--', 'color', [0.5 0.5 0.5], 'ax', ax6);
end