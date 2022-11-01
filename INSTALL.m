function INSTALL()
    % Function to install Combustion Toolbox as a package that can be used
    % using the desktop environment (plain code) or the GUI.
    % Install the MATLAB toolbox from the *.mlappinstall file

    % Get the directory path
    dir_name = fileparts( mfilename('fullpath') );
    % Find *.mlappinstall file in the installer folder
    toolbox_filename = dir(fullfile(dir_name, 'installer', '*.mltbx'));
    % Install the toolbox
    toolbox = matlab.addons.toolbox.installToolbox(fullfile(toolbox_filename.folder, toolbox_filename.name));
    % Print installed addons to the screen to verify the installation
    matlab.addons.installedAddons()
end