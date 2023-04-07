function INSTALL(varargin)
    % This function installs the Combustion Toolbox repository from local
    % files. It adds all subfolders to the MATLAB path and installs the
    % Combustion Toolbox app.
    %
    % Optional Args:
    %     * action (char): 'install' or 'uninstall' (default: 'install')
    %     * type (char): 'path', 'GUI', or 'all' (default: 'all')
    %
    % Examples:
    %     * INSTALL();                  % Installs the Combustion Toolbox
    %     * INSTALL('uninstall');       % Uninstalls the Combustion Toolbox
    %     * INSTALL('install', 'path'); % Installs the Combustion Toolbox (only MATLAB path)
    %     * INSTALL('install', 'GUI');  % Installs the Combustion Toolbox (only GUI)
    %     * INSTALL('install', 'all');  % Installs the Combustion Toolbox
    %
    % Notes:
    %     The code is available in:
    %     * GitHub - https://github.com/AlbertoCuadra/combustion_toolbox
    %     * MATLAB File Exchange - https://in.mathworks.com/matlabcentral/fileexchange/101088-combustion-toolbox
    %     * Zenodo - https://doi.org/10.5281/zenodo.5554911
    %
    % Website: https://combustion-toolbox-website.readthedocs.io/ 
    %          or type "run website_CT"
    %
    % @author: Alberto Cuadra Lara
    %          PhD Candidate - Group Fluid Mechanics
    %          Universidad Carlos III de Madrid
    
    % Default
    action = 'install';
    type = 'all';
    
    % Definitions
    name = 'Combustion Toolbox';
    app_name = 'combustion_toolbox_app';
    dir_code = get_dir_code();
    dir_app = fullfile(dir_code, 'installer', [app_name, '.mlappinstall']);

    % Unpack
    if nargin
        action = varargin{1};
    end

    if nargin > 1
        type = varargin{2};
    end
    
    % Get type of installation/uninstallation
    [FLAG_PATH, FLAG_GUI] = get_type(type);
    
    % Install/Uninstall Combustion Toolbox
    action_code(action);
    
    % NESTED FUNCTIONS
    function action_code(action)
        % Install or uninstall the Combustion Toolbox
        %
        % Args:
        %     action (char): 'install' or 'uninstall'
    
        % Definitions
        switch lower(action)
            case 'install'
                f_path = @addpath;
                f_app = @matlab.apputil.install;
                message = 'Installing';
                
            case 'uninstall'
                f_path = @rmpath;
                f_app = @matlab.apputil.uninstall;
                app_info = matlab.apputil.getInstalledAppInfo;
                FLAG_ID = contains(struct2table(app_info).id, app_name);
    
                if ~sum(FLAG_ID)
                    dir_app = [];
                else
                    dir_app =  app_info(FLAG_ID).id;
                end
    
                message = 'Uninstalling';

            otherwise
                error('Invalid action specified. Please use ''install'' or ''uninstall''.');
        end
        
        % Install/Uninstall path
        action_path(f_path, message);
        
        % Install/Uninstall the Combustion Toolbox app
        action_app(f_app, message);
    
    end
    
    function action_path(f_path, message)
        % Install/Uninstall path
        %
        % Args:
        %     f_path (function): Function to add or remove the path
        %     message (char): Message to display

        if ~FLAG_PATH
            return
        end

        % Add the code directory to the MATLAB path
        fprintf('%s %s path... ', message, name);
        f_path(dir_code);
        
        % Get subfolders (except '.git' and '.github')
        d = dir(dir_code);
        subfolders = d([d(:).isdir]);
        subfolders = subfolders(~ismember({subfolders(:).name}, {'.', '..', '.git', '.github'}));
        
        % Add all subfolders to the MATLAB path
        for i = length(subfolders):-1:1
            dir_subfolders = fullfile(subfolders(i).folder, subfolders(i).name);
            genpath_subfolders = genpath(dir_subfolders);
            f_path(genpath_subfolders);
        end
        
        fprintf('OK!\n')
    end

    function action_app(f_app, message)
        % Install/Uninstall the Combustion Toolbox app
        %
        % Args:
        %     f_app (function): Function to install or uninstall the app
        %     message (char): Message to display
        if ~FLAG_GUI
            return
        end
            
        if isempty(dir_app)
            fprintf('%s app is not currently installed.\n', name);
            return
        end

        fprintf('%s %s app...  ', message, name);
        f_app(dir_app);
        fprintf('OK!\n')
        
    end

end

% SUB-PASS FUNCTIONS
function dir_code = get_dir_code()
    % Get the directory where the code is located
    %
    % Returns:
    %     dir_code (char): Directory where the code is located
    
    % Check if user is running an executable in standalone mode
    if isdeployed
        [~, result] = system('set PATH');
        dir_code = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
        return
    end
    
    % Get current folder if the user is running an m-file from the
    % MATLAB integrated development environment (regular MATLAB)
    dir_code = pwd;
end

function [FLAG_PATH, FLAG_GUI] = get_type(type)
    % Get the type of installation/uninstallation
    %
    % Args:
    %     type (char): 'path', 'GUI', or 'all'
    %
    % Returns:
    %     Tuple containing
    %
    %     * FLAG_PATH (bool): True if the path is installed/uninstalled
    %     * FLAG_GUI (bool): True if the GUI is installed/uninstalled

    switch lower(type)
        case 'all'
            FLAG_PATH = true;
            FLAG_GUI = true;
        case 'path'
            FLAG_PATH = true;
            FLAG_GUI = false;
        case 'gui'
            FLAG_PATH = false;
            FLAG_GUI = true;
        otherwise
            error('Invalid type specified. Please use ''path'', ''gui'', or ''all''.');
    end

end