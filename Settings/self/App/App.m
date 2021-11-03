function app = App(varargin)
    % Function that create a struct called app with all the data needed
    [app, LS, FLAG_FAST] = initialize(varargin);
    app.E = Elements();
    app.S = Species();
    app.C = Constants();
    app.Misc = Miscellaneous();
    app.PD = ProblemDescription();
    app.PS = ProblemSolution();
    app.TN = TunningProperties();
    app = constructor(app, LS, FLAG_FAST);
    if ~nargin || ~isa(varargin{1,1}, 'combustion_toolbox_app') || ~isa(varargin{1,1}, 'combustion_toolbox_app_develop') || ~isa(varargin{1,1}, 'combustion_toolbox_app_original') || ~strcmpi(varargin{1,1}, 'fast')
        app = Initialize(app);
    end
end

function [app, LS, FLAG_FAST] = initialize(varargin)
    varargin = varargin{1,1};
    nargin = numel(varargin);
    app = struct();
    LS = [];
    if nargin
        if ~strcmpi(varargin{1,1}, 'fast')
            FLAG_FAST = false;
        else
            FLAG_FAST = true;
            app.DB_master = varargin{1,2};
            app.DB = varargin{1,3};
            return
        end
        if isa(varargin{1,1}, 'combustion_toolbox_app') || isa(varargin{1,1}, 'combustion_toolbox_app_develop') || isa(varargin{1,1}, 'combustion_toolbox_app_original')
            app = varargin{1,1};
            if nargin == 2
                LS = varargin{1,2};
            end
        else
            LS = varargin{1,1};
        end
    end
end

function app = constructor(app, LS, FLAG_FAST)
    % FLAG_GUI
    app = check_GUI(app);
    % Set Database
    reducedDB = false;
    app = set_DB(app, reducedDB, FLAG_FAST);
    % Set ListSpecies
    if ~isempty(LS)
        app = ListSpecies(app, LS);
    else
        app = ListSpecies(app);
    end
    % Set Contained elements
    app = ContainedElements(app);
    % Timer
    app.Misc.timer_0 = tic;
end

function app = check_GUI(app)
    if isa(app, 'combustion_toolbox_app') || isa(app, 'combustion_toolbox_app_develop') || isa(app, 'combustion_toolbox_app_original')
        app.Misc.FLAG_GUI = true;
    end
end