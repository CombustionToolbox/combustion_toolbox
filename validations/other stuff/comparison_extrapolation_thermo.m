function comparison_extrapolation_thermo(varargin)
    % Comparison of the thermodynamic functions with the higher order terms of
    % the polynomials fits and with a linearized extrapolation
    %
    % Optional Args:
    %     * species (cell): List of species
    %     * properties (cell): List of thermodynamic properties [cp, cv, h, s, gamma, g]
    %     * T (float): Temperature range
    %
    % Example:
    %     comparison_extrapolation_thermo('species', {'O2', 'N2', 'eminus', 'CO2'}, 'properties', {'cp', 'h', 's'}, 'T', 0:50:40000)
    
    % Default values
    species = {'O2', 'N2', 'eminus', 'CO2'};
    properties = {'cp', 'h', 's'};
    T = 0:50:40000;

    % Unpack
    for i = 1:2:nargin
        switch lower(varargin{i})
            case {'species', 'ls'}
                species = varargin{i + 1};
            case {'properties', 'prop', 'thermo'}
                properties = varargin{i + 1};
            case {'temperature', 't', 'temp'}
                T = varargin{i + 1};
        end

    end
    
    % Definitions
    NS = length(species);
    NP = length(properties);
    NT = length(T);
    
    % Initialization
    self = App();
    
    % Compute thermodynamic polynomials
    [values_NASA, values_CT, y_labelname] = get_thermo(self, species, properties, T, NS, NP, NT);
    
    % Get FLAG if the calculations imply extrapolation
    [FLAG_EXTRAPOLATION_PRE, FLAG_EXTRAPOLATION_POST] = get_FLAG_EXTRAPOLATION(self, species, T, NS, NT);
    ALL_FLAGS = FLAG_EXTRAPOLATION_POST & FLAG_EXTRAPOLATION_PRE;
    
    % Define color palette
    color_palette = brewermap(NS, self.Misc.config.colorpalette);
    
    % Define normalized size of the figure
    self.Misc.config.innerposition = [0.1, 0.2, 0.8, 0.7];
    self.Misc.config.outerposition = [0.1, 0.2, 0.8, 0.7];
    
    % Plot results
    set_figure(self.Misc.config);
    
    tiledlayout(1, NP);
    
    for j = 1:NP
        ax = nexttile;
        self.Misc.config.xscale = 'log';
        self.Misc.config.yscale = 'log';
        ax = set_figure(ax, self.Misc.config);

        for i = 1:NS
            temp = reshape(values_CT(i, j, :), 1, NT);
            ax = plot_figure('T', T(~ALL_FLAGS(i, :)), y_labelname{j}, temp(~ALL_FLAGS(i, :)), 'color', color_palette(i, :), 'linestyle', '-', 'ax', ax);
            ax = plot_figure('T', T(ALL_FLAGS(i, :)), y_labelname{j}, temp(ALL_FLAGS(i, :)), 'color', color_palette(i, :), 'linestyle', ':', 'ax', ax);
            plot_figure('T', T, y_labelname{j}, reshape(values_NASA(i, j, :), 1, NT), 'color', color_palette(i, :), 'linestyle', '--', 'ax', ax);
        end
    
    end

end

% SUB-PASS FUNCTIONS
function [values_NASA, values_CT, y_labelname] = get_thermo(self, species, properties, T, NS, NP, NT)
    % Compute thermodynamic polynomials
    for i = NS:-1:1
        for j = NP:-1:1
            for k = NT:-1:1
                % Get functions handles
                [funname_NASA, funname_CT, y_labelname{i}] = set_inputs_thermo_validations(properties{j});
    
                % Compute NASA
                values_NASA(i, j, k) = funname_NASA(species{i}, T(k), self.DB);
        
                % Compute CT
                values_CT(i, j, k) = funname_CT(species{i}, T(k), self.DB);
            end
    
        end
    
    end
end

function [FLAG_EXTRAPOLATION_PRE, FLAG_EXTRAPOLATION_POST] = get_FLAG_EXTRAPOLATION(self, species, T, NS, NT)
    % Get FLAG if the calculations imply extrapolation 
    FLAG_EXTRAPOLATION_PRE = false(NS, NT);
    FLAG_EXTRAPOLATION_POST = false(NS, NT);
    for i = NS:-1:1
        FLAG_EXTRAPOLATION_PRE(i, :) = T < self.DB.(species{i}).T(1);
        FLAG_EXTRAPOLATION_POST(i, :) = T > self.DB.(species{i}).T(end);
    end

end