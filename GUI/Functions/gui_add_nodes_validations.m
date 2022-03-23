function gui_add_nodes_validations(obj, code_validation_name)
    % Add nodes with the name of the validations routines
    try
        mfiledir = fileparts(fileparts(fileparts(mfilename('fullpath'))));
        filenames = dir(fullfile(mfiledir, 'Validations', code_validation_name, '*.m'));
        for i = 1:length(filenames)
            uitreenode(obj.(code_validation_name), 'Text', filenames(i).name);
        end
    catch ME
        % Print error
        fprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
        ME.stack(1).name, ME.stack(1).line, ME.message);
    end
end