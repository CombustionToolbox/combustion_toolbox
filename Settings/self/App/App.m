function app = App(varargin)
    [app, LS] = initialize(varargin);
    app.E = Elements();
    app.S = Species();
    app.C = Constants();
    app.Misc = Miscellaneous();
    app.PD = ProblemDescription();
    app.PS = ProblemSolution();
    app.TN = TunningProperties();
    app = constructor(app, LS);
    if ~isa(varargin{1,1}, 'combustion_toolbox_app') || ~isa(varargin{1,1}, 'combustion_toolbox_app_exported') || ~isa(varargin{1,1}, 'combustion_toolbox_app_original')
        app = Initialize(app);
    end
end

function [app, LS] = initialize(varargin)
    varargin = varargin{1,1};
    nargin = numel(varargin);
    app = struct();
    LS = [];
    if nargin
        if isa(varargin{1,1}, 'combustion_toolbox_app') || isa(varargin{1,1}, 'combustion_toolbox_app_exported') || isa(varargin{1,1}, 'combustion_toolbox_app_original')
            app = varargin{1,1};
            if nargin == 2
                LS = varargin{1,2};
            end
        else
            LS = varargin{1,1};
        end
    end
end

function app = constructor(app, LS)
    % FLAG_GUI
    app = check_GUI(app);
    % Set Database
    reducedDB = false;
    app = set_Database(app, reducedDB);
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
    if isa(app, 'combustion_toolbox_app') || isa(app, 'combustion_toolbox_app_exported') || isa(app, 'combustion_toolbox_app_original')
        app.Misc.FLAG_GUI = true;
    end
end