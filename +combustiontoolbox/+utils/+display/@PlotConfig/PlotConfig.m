classdef PlotConfig < handle
    % Class to store plot configuration data
    %
    % Attributes:
    %     position (float): Default figure position [pixels]
    %     innerposition (float): Set figure inner position [normalized]
    %     outerposition (float): Set figure outer position [normalized]
    %     linestyle (char): Set line style for plots
    %     symbolstyle (char): Set symbol style for plots
    %     linewidth (float): Set line width for plots
    %     fontsize (float): Set fontsize
    %     colorpalette (char): Set Color palette (see brewermap function for more options)
    %     colorpaletteLenght (float): Set Maximum number of colors to use in the color palette
    %     box (char): Display axes outline
    %     grid (char): Display or hide axes grid lines
    %     hold (char): Retain current plot when adding new plots
    %     axis_x (char): Set x-axis limits
    %     axis_y (char): Set y-axis limits
    %     xscale (char): Set x-axis scale (linear or log)
    %     yscale (char): Set y-axis scale (linear or log)
    %     xdir (char): Set x-axis direction (normal or reverse)
    %     ydir (char): Set y-axis direction (normal or reverse)
    %     title (char): Set title
    %     label_type (char): Set label with variable (short), name (medium), or name and variable (long)
    %     labelx (char): Set x label
    %     labely (char): Set y label
    %     legend_name (char): Set legend labels
    %     legend_location (char): Set legend location
    %     colorline (float): Default colorline
    %     colorlines (float): Default colorlines
    %     blue (float): Default colorline
    %     gray (float): Default colorline
    %     red (float): Default colorline
    %     orange (float): Default colorline
    %     brown (float): Default colorline
    %     brown2 (float): Default colorline
    %     id_polar1 (float): Axes id for pressure-deflection polar diagrams
    %     id_polar2 (float): Axes id for wave-deflection polar diagrams
    %     id_polar3 (float): Axes id for velocity polar diagrams
    
    properties
        innerposition = [0.05 0.05 0.9 0.9]   % Set figure inner position [normalized]
        outerposition = [0.05 0.05 0.9 0.9]   % Set figure outer position [normalized]
        linestyle = '-'                       % Set line style for plots
        symbolstyle = 'o'                     % Set symbol style for plots
        linewidth = 1.8                       % Set line width for plots
        fontsize = 22                         % Set fontsize
        colorpalette = 'Seaborn'              % Set Color palette (see brewermap function for more options)
        colorpaletteLenght = 11               % Set Maximum number of colors to use in the color palette
        box = 'off'                           % Display axes outline
        grid = 'off'                          % Display or hide axes grid lines
        hold = 'on'                           % Retain current plot when adding new plots
        axis_x = 'tight'                      % Set x-axis limits
        axis_y = 'auto'                       % Set y-axis limits
        xscale = 'linear'                     % Set x-axis scale (linear or log)
        yscale = 'linear'                     % Set y-axis scale (linear or log)
        xdir = 'normal'                       % Set x-axis direction (normal or reverse)
        ydir = 'normal'                       % Set y-axis direction (normal or reverse)
        title = []                            % Set title
        label_type = 'short'                  % Set label with variable (short), name (medium), or name and variable (long)
        labelx = []                           % Set x label
        labely = []                           % Set y label
        legend_name = []                      % Set legend labels
        legend_location = 'northeastoutside'  % Set legend location
        colorline = [44, 137, 160]/255        % Default colorline
        colorlines = [135, 205, 222;...       % Default colorlines
                      95, 188, 211;...        % Default colorlines
                      44, 137, 160;...        % Default colorlines
                      22,  68,  80]/255       % Default colorlines
        blue = [0.3725, 0.7373, 0.8275]       % Default colorline
        gray = [0.50, 0.50, 0.50]             % Default colorline
        red = [0.64,0.08,0.18]                % Default colorline
        orange = [212, 85, 0]/255             % Default colorline
        brown = [200, 190, 183]/255           % Default colorline
        brown2 = [72, 55, 55]/255             % Default colorline
        id_polar1 = 1001                      % Axes id for pressure-deflection polar diagrams
        id_polar2 = 1002                      % Axes id for wave-deflection polar diagrams
        id_polar3 = 1003                      % Axes id for velocity polar diagrams
        displaySpecies                        % Display species
        mintolDisplay = 1e-6                  % Minimum tolerance to display species
        plotProperties = {'T', 'p', 'rho', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound'} % Plot properties
        plotPropertiesBasis = {[], [], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], []} % Plot properties basis
    end

    properties (Dependent)
        position          % Figure position [pixels]
        numPlotProperties % Number of properties to plot
    end

    properties (Access = private)
        position_
    end

    methods
        function value = get.numPlotProperties(obj)
            value = length(obj.plotProperties);
        end

        function value = get.position(obj)
            % Get figure position
            if ~isempty(obj.position)
                value = combustiontoolbox.utils.display.getMonitorPositions(2);
                obj.position_ = value;
                return
            end

            value = obj.position_;
        end

        function set.position(obj, value)
            % Set figure position
            if ~isempty(value)
                obj.position_ = value;
            end
            
        end

    end

    methods
        
        function obj = PlotConfig(varargin)
            % Class constructor

            % Set properties
            if nargin > 0
                obj = set(obj, varargin{:});
            end

        end

        function obj = set(obj, property, value, varargin)
            % Set properties of the PlotConfig object
            %
            % Args:
            %     obj (PlotConfig): PlotConfig object
            %     property (char): Property name
            %     value (float): Property value
            %
            % Optional Args:
            %     * property (char): Property name
            %     * value (float): Property value
            %
            % Returns:
            %     obj (PlotConfig): PlotConfig object with updated properties
            %
            % Examples:
            %     * set(PlotConfig(), 'fontsize', 18)
            %     * set(PlotConfig(), 'fontsize', 18, 'linewidth', 2)
            
            varargin = [{property, value}, varargin{:}];

            for i = 1:2:length(varargin)
                % Assert that the property exists
                assert(isprop(obj, varargin{i}), 'Property not found');

                % Set property
                obj.(varargin{i}) = varargin{i + 1};
            end

        end

    end

end