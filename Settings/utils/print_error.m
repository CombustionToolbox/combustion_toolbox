function error_message = print_error(ME, varargin)
    % Print message error
    %
    % Args:
    %     ME (obj): MException object that allows to identify the error
    %
    % Optional Args:
    %     message_solution (str): Message solution
    %
    % Reutnrs:
    %     error_message (str): Message error

    % Default
    type = 'Error';
    message_solution = 'None.';
    % Unpack
    for i = 1:2:nargin - 1

        switch lower(varargin{i})
            case {'solution', 'sol', 'action'}
                message_solution = varargin{i + 1};
            case {'type'}
                type = varargin{i + 1};
        end

    end
    
    % Set error message
    try
        error_message = sprintf('%s in function %s() at line %d.\n%s Message: %s', ...
            type, ME.stack(2).name, ME.stack(2).line, type, ME.message);
    catch
        error_message = sprintf('%s in function %s() at line %d.\n%s Message: %s', ...
            type, ME.stack(1).name, ME.stack(1).line, type, ME.message);
    end
    
    % Set function
    switch lower(type)
        case 'error'
            f_print = @error;
        case 'warning'
            f_print = @warning;
        otherwise
            f_print = @fprintf;
    end
    % Set error solution
    error_solution = sprintf('%s Solution: %s', type, message_solution);
    % Print message
    f_print('%s\n%s\n', error_message, error_solution);
end
