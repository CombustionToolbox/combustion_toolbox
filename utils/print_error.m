function error_message = print_error(ME, varargin)
    % Print message error
    %
    % Args:
    %     ME (object): MException object that allows to identify the error
    %
    % Optional Name-Value Pairs Args:
    %     * type (char): Type of message (error, warning, or other)
    %     * message_solution (char): Message solution
    %
    % Returns:
    %     error_message (char): Message error
    %
    % Examples:
    %     * error_message = print_error(ME, 'Type', 'Warning')
    %     * error_message = print_error(ME, 'Type', 'Warning', 'Solution', 'Returning an empty index value.')

    % Default
    type = 'Error';
    message_solution = 'None.';

    % Unpack
    for i = 1:2:nargin - 1

        switch lower(varargin{i})
            case {'type'}
                type = varargin{i + 1};
            case {'solution', 'sol', 'action'}
                message_solution = varargin{i + 1};
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
