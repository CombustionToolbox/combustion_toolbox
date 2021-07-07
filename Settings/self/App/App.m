function app = App(varargin)
    [app, minors] = initialize(varargin);
    app.E = Elements();
    app.S = Species();
    app.M = MinorsProducts();
    app.C = Constants();
    app.Misc = Miscelaneous();
    app.PD = ProblemDescription();
    app.PS = ProblemSolution();
    app.TN = TunningProperties();
    app = constructor(app, minors);
end

function [app, minors] = initialize(varargin)
    varargin = varargin{1,1};
    nargin = numel(varargin);
    app = struct();
    minors = [];
    if nargin
        if isa(varargin{1,1}, 'combustion_toolbox_app')
            app = varargin{1,1};
            if nargin == 2
                minors = varargin{1,2};
            end
        else
            minors = varargin{1,1};
        end
    end
end

function app = constructor(app, minors)
    % Timer
    app.Misc.timer_0 = tic;
    % FLAG_GUI
    app = check_GUI(app);
    % Set Database
    reducedDB = true;
    app = set_Database(app, reducedDB);
    % Contained elements
    app = ContainedElements(app);
    % Definition MinorsProducts
    if ~isempty(minors)
        app = MinorsProducts_self(app, minors);
    else
        app = MinorsProducts_self(app);
    end
end

function app = check_GUI(app)
    if isa(app, 'combustion_toolbox_app')
        app.Misc.FLAG_GUI = true;
    end
end