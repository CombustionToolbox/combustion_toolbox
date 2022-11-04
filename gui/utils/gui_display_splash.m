function splash_obj = gui_display_splash(varargin)
    % Display splash using SplashScreen [1]
    %
    % [1] Ben Tordoff (2022). SplashScreen (https://www.mathworks.com/matlabcentral/fileexchange/30508-splashscreen), MATLAB Central File Exchange.
    %
    % Args:
    %    app (object): App object
    %
    % Returns:
    %    splash_obj (object): Splash object
    
    % Default values
    color_splash = [0.5098, 0.6039, 0.6745];
    time_pause = [];
    % Unpack
    for i = 1:2:nargin
        switch lower(varargin{i})
            case 'app'
                color_splash = varargin{i+1}.color_splash; % [r g b] from 0 to 1
            case 'color'
                color_splash = varargin{i+1}; % [r g b] from 0 to 1
            case 'pause'
                time_pause = varargin{i+1}; % [seconds]
        end
    end
    % Get Combustion Toolbox version and published date
    [release, date] = get_combustion_toolbox_version();
    % Generate splash
    splash_obj = SplashScreen('Combustion Toolbox', fullfile('gui', 'assets', 'splash', 'splash.png'));               
    splash_obj.addText(110, 100, [release,' (', date,')'], 'FontSize', 18, 'Color', color_splash, 'FontName', 'Arial', 'Shadow', 'off');
    % Optional: auto-delete after x seconds
    if ~isempty(time_pause)
        pause(time_pause);
        delete(splash_obj);
        splash_obj = [];
    end
end