function INSTALL(varargin)
    % This function installs the Combustion Toolbox repository from local
    % files. It adds all subfolders to the MATLAB path and installs the
    % Combustion Toolbox app.
    %
    % Optional Args:
    %     * action (char): 'install' or 'uninstall' (default: 'install')
    %     * type (char): 'path', 'GUI', or 'all' (default: 'all')
    %     * FLAG_HOME (bool): true or false (default: false)
    %
    % Examples:
    %     * INSTALL();                       % Installs the Combustion Toolbox
    %     * INSTALL('uninstall');            % Uninstalls the Combustion Toolbox
    %     * INSTALL('install', 'path');      % Installs the Combustion Toolbox (only MATLAB path)
    %     * INSTALL('install', 'GUI');       % Installs the Combustion Toolbox (only GUI)
    %     * INSTALL('install', 'all');       % Installs the Combustion Toolbox
    %     * INSTALL('install', 'all', true); % Installs the Combustion Toolbox with FLAG_HOME set to true
    %
    % Notes:
    %     The code is available in:
    %     * GitHub - https://github.com/CombustionToolbox/combustion_toolbox
    %     * MATLAB File Exchange - https://in.mathworks.com/matlabcentral/fileexchange/101088-combustion-toolbox
    %     * Zenodo - https://doi.org/10.5281/zenodo.5554911
    %
    % Website: https://combustion-toolbox-website.readthedocs.io/ 
    %          or type "run website_CT"
    %
    % @author: Alberto Cuadra Lara
    %          Postdoctoral researcher - Group Fluid Mechanics
    %          Universidad Carlos III de Madrid
    
    % Default
    action = 'install';
    type = 'all';
    FLAG_HOME = false;

    % Unpack
    if nargin
        action = varargin{1};
    end

    if nargin > 1
        type = varargin{2};
    end

    if nargin > 2
        FLAG_HOME = varargin{3};
    end
    
    % Get package destination
    packageDst = getPath(FLAG_HOME);
    
    % Install Combustion Toolbox package and app
    installPackage(action, type, packageDst);
end

% SUB-PASS FUNCTIONS
function packageDst = getPath(FLAG_HOME)
    % Get path package destination

    % Get path package destination based on current location
    if ~FLAG_HOME
        packageDst = pwd();
        return
    end
    
    % Get path package destination based on home operating system location
    switch computer('arch')
        case {'mac', 'maci64', 'maca64', 'maci'}
            % macOS
            packageDst = '/Applications/CombustionToolbox/matlab';
        case 'linux64'
            % Linux
            packageDst = '/usr/local/CombustionToolbox/matlab';
        case 'win64'
            % Windows
            packageDst = 'C:\Program Files\CombustionToolbox\matlab';
        otherwise
            error('Unsupported operating system');
    end

end

function bashScriptFile = generateBash(packageDst, action)
    % Generate bash script

    % Definitions
    packageSrc = pwd();
    startupFile = fullfile(getenv('HOME'), 'Documents', 'MATLAB', 'startup.m');
    addPathLine = sprintf('addpath(genpath(''%s'')) %%%% added by Combustion toolbox installer', packageDst);
    
    % Define the content of the Bash script
    if strcmp(action, 'install')
        bashScriptContent = sprintf([
            '#!/bin/bash\n' ...
            '\n' ...
            'PACKAGE_SRC="%s"\n' ...
            'PACKAGE_DST="%s"\n' ...
            'STARTUP_FILE="%s"\n' ...
            'ADDPATH_LINE="%s"\n' ...
            '\n' ...
            'printf "Copying Combustion Toolbox destination directory... "\n' ...
            'rsync -aq --exclude=".git" "$PACKAGE_SRC/" "$PACKAGE_DST/"\n' ...
            'printf "OK!\\n"\n' ...
            '\n' ...
            'printf "Adding Combustion Toolbox destination directory to startup.m... "\n' ...
            'if grep -Fxq "$ADDPATH_LINE" "$STARTUP_FILE"; then\n' ...
            '    printf "OK!\\n"\n' ...
            'else\n' ...
            '    echo "$ADDPATH_LINE" >> "$STARTUP_FILE"\n' ...
            '    printf "OK!\\n"\n' ...
            'fi\n'], ...
        packageSrc, packageDst, startupFile, addPathLine);
    elseif strcmp(action, 'uninstall')
        bashScriptContent = sprintf([
            '#!/bin/bash\n' ...
            '\n' ...
            'STARTUP_FILE="%s"\n' ...
            '\n' ...
            'printf "Removing Combustion Toolbox data from startup.m... "\n' ...
            'sed -i '''' "/added by Combustion toolbox installer/d" "$STARTUP_FILE"\n' ...
            'printf "OK!\\n"\n'], ...
        startupFile);
    else
        error('Invalid action specified. Please use ''install'' or ''uninstall''.');
    end
    
    % Write the Bash script content to a file
    bashScriptFile = 'install.sh';
    fid = fopen(bashScriptFile, 'wt');
    fprintf(fid, bashScriptContent);
    fclose(fid);

    % Make the script executable
    system(['chmod +x ' bashScriptFile]);
end

function installPackage(action, type, packageDst)
    % This function installs the Combustion Toolbox repository from local
    % files. It adds all subfolders to the MATLAB path and installs the
    % Combustion Toolbox app.
    %
    % Args:
    %     action (char): 'install' or 'uninstall'
    %     type (char): 'path', 'GUI', or 'all'
    %     packageDst (char): package destination path
    %
    % Examples:
    %     * installPackage();                  % Installs the Combustion Toolbox
    %     * installPackage('uninstall');       % Uninstalls the Combustion Toolbox
    %     * installPackage('install', 'path'); % Installs the Combustion Toolbox (only MATLAB path)
    %     * installPackage('install', 'GUI');  % Installs the Combustion Toolbox (only GUI)
    %     * installPackage('install', 'all');  % Installs the Combustion Toolbox
    %
    % Notes:
    %     The code is available in:
    %     * GitHub - https://github.com/CombustionToolbox/combustion_toolbox
    %     * MATLAB File Exchange - https://in.mathworks.com/matlabcentral/fileexchange/101088-combustion-toolbox
    %     * Zenodo - https://doi.org/10.5281/zenodo.5554911
    %
    % Website: https://combustion-toolbox-website.readthedocs.io/ 
    %          or type "run website_CT"
    %
    % @author: Alberto Cuadra Lara
    %          Postdoctoral researcher - Group Fluid Mechanics
    %          Universidad Carlos III de Madrid
        
    % Definitions
    name = 'Combustion Toolbox';
    app_name = 'combustion_toolbox_app';
    dir_code = getDirCode();
    dir_database = fullfile(dir_code, 'databases');
    dir_app = fullfile(dir_code, 'installer', [app_name, '.mlappinstall']);
    
    % Get type of installation/uninstallation
    [FLAG_PATH, FLAG_GUI] = getType(type);
    
    % Add/remove Combustion Toolbox in MATLAB startup.m file
    try
        % Generate bash script
        bashScriptFile = generateBash(packageDst, action);
        
        % Execute bash script
        system(['./' bashScriptFile]);
    
        % Remove the bash script
        delete(bashScriptFile);
    catch
        fprintf('Error including Combustion Toolbox path in MATLAB startup.m file');
    end

    % Move local path to package destination
    cd(packageDst);

    % Install/Uninstall Combustion Toolbox
    actionCode(action);
    
    % NESTED FUNCTIONS
    function actionCode(action)
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

                if isempty(app_info)
                    return
                end
                
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
        actionPath(f_path, message);

        % Install/Uninstall the Combustion Toolbox app
        actionApp(f_app, message);
    end
    
    function actionDatabases(action, path)
        % Install/Uninstall databases
        %
        % Args:
        %     action (char): 'install' or 'uninstall'

        % Import packages
        import combustiontoolbox.databases.*

        if strcmpi(action, 'uninstall')
            return
        end

        % Install databases
        databases = {NasaDatabase()};

        for i = 1:length(databases)
            databases{i}.save('path', path);
        end

    end

    function actionPath(f_path, message)
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

            % Check if the subfolder is a package folder
            if startsWith(subfolders(i).name, '+')
                f_path(subfolders(i).folder); % Add parent directory
                continue
            end

            genpath_subfolders = genpath(dir_subfolders);
            f_path(genpath_subfolders);
        end
        
        % Save the path permanently
        if strcmpi(action, 'install')
            savepath;

            % Install databases
            actionDatabases(action, dir_database);
        end

        fprintf('OK!\n')
    end

    function actionApp(f_app, message)
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
        
        % Get dir app
        app_info = matlab.apputil.getInstalledAppInfo;
        FLAG_ID = contains(struct2table(app_info).id, app_name);
        dir_app =  app_info(FLAG_ID).location;
        dir_database_app = fullfile(dir_app, 'databases');

        % Make directory databases
        if ~exist(dir_database_app, 'dir')
            mkdir(dir_database_app);
            % Add directory to the MATLAB path
            addpath(dir_database_app);
            % Save the path permanently
            savepath;
        end
        
        % Install databases
        actionDatabases(action, dir_database_app);
    end

end

function dir_code = getDirCode()
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

function [FLAG_PATH, FLAG_GUI] = getType(type)
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