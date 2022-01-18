function Combustion_Toolbox_setPath()
    % Set directories to MATLAB search path
    
    % Get current path
    executableFolder = GetExecutableFolder();
    % Generate a path that includes myfolder and all folders below it
    p = genpath(executableFolder);
    % Add the folder and its subfolders to the search path
    addpath(p);
end

% SUB-PASS FUNCTIONS
function [executableFolder] = GetExecutableFolder()
    % Returns the folder where the compiled executable actually resides
    try
        if isdeployed
            % User is running an executable in standalone mode.
            [status, result] = system('set PATH');
            executableFolder = char(regexpi(result, 'Path=(.*?);', 'tokens', 'once'));
            % 			fprintf(1, '\nIn function GetExecutableFolder(), currentWorkingDirectory = %s\n', executableFolder);
        else
            % User is running an m-file from the MATLAB integrated development environment (regular MATLAB).
            executableFolder = pwd;
        end
    catch ME
        errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
            ME.stack(1).name, ME.stack(1).line, ME.message);
        uiwait(warndlg(errorMessage));
    end
end