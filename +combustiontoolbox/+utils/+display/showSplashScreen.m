function splashScreen = showSplashScreen(varargin)
    % Display the Combustion Toolbox splash screen
    %
    % Optional Name-Value Args:
    %     app (object): Combustion Toolbox app object
    %     color (double): Normalized RGB text color
    %     pause (double): Time to pause before deleting the splash [s]
    %
    % Returns:
    %     splashScreen (SplashScreen): Startup splash screen object

    color = [0.5098, 0.6039, 0.6745];
    pauseTime = [];

    for i = 1:2:nargin
        name = lower(char(varargin{i}));
        value = varargin{i + 1};

        switch name
            case 'app'
                color = appSplashColor(value, color);
            case 'color'
                color = value;
            case 'pause'
                pauseTime = value;
        end
    end

    constants = combustiontoolbox.common.Constants;
    imagePath = splashImagePath();
    splashScreen = combustiontoolbox.utils.extensions.SplashScreen( ...
        'Combustion Toolbox', imagePath);
    splashScreen.addText(110, 100, [constants.release, ' (', constants.date, ')'], ...
        'FontSize', 18, 'Color', color, 'FontName', 'Arial', 'Shadow', 'off');

    if ~isempty(pauseTime)
        pause(pauseTime);
        delete(splashScreen);
        splashScreen = [];
    end
end

function color = appSplashColor(app, defaultColor)
    % Return splash color from the app when available
    color = defaultColor;

    if isobject(app) && isprop(app, 'splashColor') && ~isempty(app.splashColor)
        color = app.splashColor;
    end
end

function path = splashImagePath()
    % Return absolute path to the splash image
    displayFolder = fileparts(mfilename('fullpath'));
    utilitiesFolder = fileparts(displayFolder);
    packageFolder = fileparts(utilitiesFolder);
    repositoryRoot = fileparts(packageFolder);
    path = fullfile(repositoryRoot, 'gui', 'assets', 'splash', 'splash.png');
end
