function splashScreen = showSplashScreen(varargin)
    % Show Combustion Toolbox splash screen
    %
    % Optional Name-Value Pairs Args:
    %     * color (float): Normalized RGB color
    %     * pause (float): Time before closing the splash screen [s]
    %
    % Returns:
    %     splashScreen (matlab.ui.Figure): Splash screen figure

    % Parse inputs
    parser = inputParser;
    parser.addParameter('color', [0.5098, 0.6039, 0.6745], @(value) isnumeric(value) && numel(value) == 3);
    parser.addParameter('pause', [], @(value) isempty(value) || (isnumeric(value) && isscalar(value) && value >= 0));
    parser.parse(varargin{:});

    color = parser.Results.color;
    timePause = parser.Results.pause;

    % Get Combustion Toolbox version and published date
    release = combustiontoolbox.common.Constants.release;
    date = combustiontoolbox.common.Constants.date;

    % Get splash path
    splashPath = getSplashPath();

    % Show splash screen
    [width, height] = getSplashSize(splashPath);
    position = getSplashPosition(width, height);
    splashScreen = uifigure(...
        'Visible', 'off', ...
        'Name', 'Combustion Toolbox', ...
        'Resize', 'off', ...
        'Color', [1, 1, 1], ...
        'Position', position);
    uiimage(splashScreen, 'ImageSource', splashPath, 'Position', [1, 1, width, height]);
    uilabel(splashScreen, ...
        'Text', [release, ' (', date, ')'], ...
        'FontName', 'Arial', ...
        'FontSize', 18, ...
        'FontColor', color, ...
        'Position', [110, height - 115, 260, 30]);
    splashScreen.Visible = 'on';
    drawnow;
    figure(splashScreen);
    drawnow;

    % Close splash after an optional delay
    if ~isempty(timePause)
        pause(timePause);

        if isvalid(splashScreen)
            delete(splashScreen);
        end

        splashScreen = [];
    end

end

function position = getSplashPosition(width, height)
    % Get centered splash figure position
    monitorPosition = combustiontoolbox.utils.display.getMonitorPositionsMATLAB();
    x0 = monitorPosition(1) + (monitorPosition(3) - width) / 2;
    y0 = monitorPosition(2) + (monitorPosition(4) - height) / 2;
    position = round([x0, y0, width, height]);
end

function splashPath = getSplashPath()
    % Get splash image path
    toolboxRoot = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    splashPath = fullfile(toolboxRoot, 'gui', 'assets', 'splash', 'splash.png');

    if exist(splashPath, 'file') == 2
        return
    end

    splashPath = fullfile('gui', 'assets', 'splash', 'splash.png');
end

function [width, height] = getSplashSize(splashPath)
    % Get splash image size
    imageInfo = imfinfo(splashPath);
    width = imageInfo.Width;
    height = imageInfo.Height;
end
