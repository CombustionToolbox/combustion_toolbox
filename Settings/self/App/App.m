function app = App(varargin)
    app.E = Elements();
    app.S = Species();
    app.M = MinorsProducts();
    app.C = Constants();
    app.Misc = Miscelaneous();
    app.PD = ProblemDescription();
    app.PS = ProblemSolution();
    app.TN = TunningProperties();
    % Constructor
    if nargin < 1
        minors = [];
    else
        minors = varargin{1,1};
    end
    app = constructor(app, minors);
end
function app = constructor(app, minors)
    % Timer
    app.Misc.timer_0 = tic;
    % Set Database
    reducedDB = true;
    app = set_Database(app, reducedDB);
    % Contained elements
    app = ContainedElements(app);
    % Definition MinorsProducts
    if minors
        app = MinorsProducts_self(app, minors);
    else
        app = MinorsProducts_self(app);
    end
end