% VALIDATION: TP_TEA
%
% Compute equilibrium composition at defined temperature and pressure.
% Reproduce the example case of TEA by Jasmina Blecic.
% URL RESULTS TEA:
% https://github.com/dzesmin/TEA/tree/master/doc/examples/quick_example/results 
%   
%
% @author: Alberto Cuadra Lara
%          PhD Candidate - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                  
% Last update April 15 2022

% Inputs
load exoplanet_TP_K2-18b.mat Pressure Temp

metallicity = 100;
Fuel = {'H', 'C', 'N', 'O'};
N_Fuel = abundances2moles(Fuel, 'abundances.txt', metallicity);

LS = {'CH4', 'CO2', 'CO', 'H2', 'H2O', 'N2', 'NH3', 'C2H2_acetylene', 'C2H6', 'HCN'};
T = linspace(Temp(1), Temp(end), 300);
p = logspace(-8, 3, 300);
Oxidizer = {};
% Tunning paramenters
mintol = 1e-32;
% Custom Plots 
DisplaySpecies = LS;
% Combustion Toolbox
results_CT = run_CT('ProblemType', 'TP',...
                    'Temperature', T,...
                    'Pressure', p, ...
                    'Species', LS,...
                    'S_Fuel', Fuel,...
                    'N_Fuel', N_Fuel,...
                    'S_Oxidizer', Oxidizer,...
                    'tolN', mintol);
problems_solved = length(results_CT.PD.range);
% Display validation (plot)
% * Molar fractions
fig1 = plot_molar_fractions_validation(results_CT, [], 'Xi', 'p', DisplaySpecies, 'mintol', mintol, 'nfrec', 3,...
        'ydir', 'reverse', 'xscale', 'log');
% Save plots
folderpath = strcat(pwd,'\Validations\Figures\');
stack_trace = dbstack;
filename = stack_trace.name;
saveas(fig1, strcat(folderpath, filename, '_molar'), 'svg');