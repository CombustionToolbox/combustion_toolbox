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
        position                              % Default figure position [pixels]
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
        displaySpecies
        mintolDisplay = 1e-6
        plotProperties = {'T', 'p', 'rho', 'h', 'e', 'g', 'cp', 's', 'gamma_s', 'sound'}
        plotPropertiesBasis = {[], [], [], 'mi', 'mi', 'mi', 'mi', 'mi', [], []}
    end

    properties (Dependent)
        numPlotProperties
    end

    methods
        function value = get.numPlotProperties(obj)
            value = length(obj.plotProperties);
        end
    end

    methods (Access = public)
        function obj = PlotConfig(varargin)
            % Class constructor
            
            % Default values
            defaultPostion = combustiontoolbox.utils.display.getMonitorPositions(2);

            % Create input parser
            ip = inputParser;
            addParameter(ip, 'position', defaultPostion, @isnumeric);
            addParameter(ip, 'innerposition', obj.innerposition, @isnumeric);
            addParameter(ip, 'outerposition', obj.outerposition, @isnumeric);
            addParameter(ip, 'linestyle', obj.linestyle, @ischar);
            addParameter(ip, 'symbolstyle', obj.symbolstyle, @ischar);
            addParameter(ip, 'linewidth', obj.linewidth, @isnumeric);
            addParameter(ip, 'fontsize', obj.fontsize, @isnumeric);
            addParameter(ip, 'colorpalette', obj.colorpalette, @ischar);
            addParameter(ip, 'colorpaletteLenght', obj.colorpaletteLenght, @isnumeric);
            addParameter(ip, 'box', obj.box, @ischar);
            addParameter(ip, 'grid', obj.grid, @ischar);
            addParameter(ip, 'hold', obj.hold, @ischar);
            addParameter(ip, 'axis_x', obj.axis_x, @ischar);
            addParameter(ip, 'axis_y', obj.axis_y, @ischar);
            addParameter(ip, 'xscale', obj.xscale, @ischar);
            addParameter(ip, 'yscale', obj.yscale, @ischar);
            addParameter(ip, 'xdir', obj.xdir, @ischar);
            addParameter(ip, 'ydir', obj.ydir, @ischar);
            addParameter(ip, 'title', obj.title, @ischar);
            addParameter(ip, 'label_type', obj.label_type, @ischar);
            addParameter(ip, 'labelx', obj.labelx, @ischar);
            addParameter(ip, 'labely', obj.labely, @ischar);
            addParameter(ip, 'legend_name', obj.legend_name, @ischar);
            addParameter(ip, 'legend_location', obj.legend_location, @ischar);
            addParameter(ip, 'colorline', obj.colorline, @isnumeric);
            addParameter(ip, 'colorlines', obj.colorlines, @isnumeric);
            addParameter(ip, 'blue', obj.blue, @isnumeric);
            addParameter(ip, 'gray', obj.gray, @isnumeric);
            addParameter(ip, 'red', obj.red, @isnumeric);
            addParameter(ip, 'orange', obj.orange, @isnumeric);
            addParameter(ip, 'brown', obj.brown, @isnumeric);
            addParameter(ip, 'brown2', obj.brown2, @isnumeric);
            addParameter(ip, 'id_polar1', obj.id_polar1, @isnumeric);
            addParameter(ip, 'id_polar2', obj.id_polar2, @isnumeric);
            addParameter(ip, 'id_polar3', obj.id_polar3, @isnumeric);
            parse(ip, varargin{:});

            % Set properties
            obj.position = ip.Results.position;
            obj.innerposition = ip.Results.innerposition;
            obj.outerposition = ip.Results.outerposition;
            obj.linestyle = ip.Results.linestyle;
            obj.symbolstyle = ip.Results.symbolstyle;
            obj.linewidth = ip.Results.linewidth;
            obj.fontsize = ip.Results.fontsize;
            obj.colorpalette = ip.Results.colorpalette;
            obj.colorpaletteLenght = ip.Results.colorpaletteLenght;
            obj.box = ip.Results.box;
            obj.grid = ip.Results.grid;
            obj.hold = ip.Results.hold;
            obj.axis_x = ip.Results.axis_x;
            obj.axis_y = ip.Results.axis_y;
            obj.xscale = ip.Results.xscale;
            obj.yscale = ip.Results.yscale;
            obj.xdir = ip.Results.xdir;
            obj.ydir = ip.Results.ydir;
            obj.title = ip.Results.title;
            obj.label_type = ip.Results.label_type;
            obj.labelx = ip.Results.labelx;
            obj.labely = ip.Results.labely;
            obj.legend_name = ip.Results.legend_name;
            obj.legend_location = ip.Results.legend_location;
            obj.colorline = ip.Results.colorline;
            obj.colorlines = ip.Results.colorlines;
            obj.blue = ip.Results.blue;
            obj.gray = ip.Results.gray;
            obj.red = ip.Results.red;
            obj.orange = ip.Results.orange;
            obj.brown = ip.Results.brown;
            obj.brown2 = ip.Results.brown2;
            obj.id_polar1 = ip.Results.id_polar1;
            obj.id_polar2 = ip.Results.id_polar2;
            obj.id_polar3 = ip.Results.id_polar3;
        end

    end

end