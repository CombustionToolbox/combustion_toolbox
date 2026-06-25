classdef AppSession < handle
    % Stores long-lived Combustion Toolbox GUI services and shared state.
    %
    % Attributes:
    %     constants (Constants): Combustion Toolbox constants
    %     database (Database): Thermodynamic database used by the GUI
    %     chemicalSystem (ChemicalSystem): Current chemical system
    %     mixture (Mixture): Current mixture object
    %     equilibriumSolver (EquilibriumSolver): CT-EQUIL solver
    %     shockSolver (ShockSolver): CT-SD shock solver
    %     detonationSolver (DetonationSolver): CT-SD detonation solver
    %     rocketSolver (RocketSolver): CT-ROCKET solver
    %     shockTurbulenceSolver (ShockTurbulenceSolver): Shock-turbulence solver
    %     plotConfig (PlotConfig): Shared plot configuration
    %     export (Export): Shared export configuration
    %
    % Examples:
    %     * session = combustiontoolbox.gui.AppSession();
    %     * session = combustiontoolbox.gui.AppSession.fromApp(app);

    properties
        constants
        database
        chemicalSystem
        mixture
        equilibriumSolver
        shockSolver
        detonationSolver
        rocketSolver
        shockTurbulenceSolver
        plotConfig
        export
    end

    methods
        function obj = AppSession(varargin)
            % AppSession constructor
            %
            % Optional Args:
            %     initialize (logical): True to create default CT services
            %
            % Returns:
            %     obj (AppSession): Initialized session object
            defaultInitialize = true;

            ip = inputParser;
            addParameter(ip, 'initialize', defaultInitialize, @islogical);
            parse(ip, varargin{:});

            if ip.Results.initialize
                initialize(obj);
            end
        end

        function initialize(obj)
            % Initialize default CT services for a standalone GUI session
            obj.constants = combustiontoolbox.common.Constants;
            obj.database = combustiontoolbox.databases.NasaDatabase();
            obj.chemicalSystem = combustiontoolbox.core.ChemicalSystem(obj.database);
            obj.mixture = combustiontoolbox.core.Mixture(obj.chemicalSystem);
            obj.plotConfig = combustiontoolbox.utils.display.PlotConfig;
            obj.equilibriumSolver = combustiontoolbox.equilibrium.EquilibriumSolver('plotConfig', obj.plotConfig);
            obj.shockSolver = combustiontoolbox.shockdetonation.ShockSolver('plotConfig', obj.plotConfig, 'equilibriumSolver', obj.equilibriumSolver);
            obj.detonationSolver = combustiontoolbox.shockdetonation.DetonationSolver('plotConfig', obj.plotConfig, 'equilibriumSolver', obj.equilibriumSolver);
            obj.rocketSolver = combustiontoolbox.rocket.RocketSolver('plotConfig', obj.plotConfig, 'equilibriumSolver', obj.equilibriumSolver);
            obj.shockTurbulenceSolver = combustiontoolbox.shockturbulence.ShockTurbulenceSolver('plotConfig', obj.plotConfig, 'equilibriumSolver', obj.equilibriumSolver, 'shockSolver', obj.shockSolver);
            obj.export = combustiontoolbox.utils.Export();
        end
    end

    methods (Static)
        function obj = fromApp(app)
            % Create a session that references services already owned by an app
            %
            % Args:
            %     app (combustion_toolbox): Main App Designer object
            %
            % Returns:
            %     obj (AppSession): Session referencing the app-owned services
            obj = combustiontoolbox.gui.AppSession('initialize', false);
            propertyNames = {'constants', 'database', 'chemicalSystem', 'mixture', ...
                'equilibriumSolver', 'shockSolver', 'detonationSolver', ...
                'rocketSolver', 'shockTurbulenceSolver', 'plotConfig', 'export'};

            for i = 1:numel(propertyNames)
                name = propertyNames{i};

                if isobject(app) && isprop(app, name)
                    obj.(name) = app.(name);
                end
            end
        end
    end
end
