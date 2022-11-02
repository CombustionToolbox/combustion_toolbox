function vector = cell2vector(value, varargin)
    % Convert values of an individual cell into a vector. If the value
    % correspond with a struct it can return as a vector the values of a
    % given fieldname.
    %
    % Args:
    %     value (cell or struct): Data of the mixture, conditions, and databases
    % Optional Args:
    %     field (str): Fieldname of the given value (struct)
    %
    % Returns:
    %     vector (any): Vector with the values of the individual cell/fieldname (struct)

    if nargin > 1
        field = varargin{1};
    end

    if nargin > 2
        error('Error in @cell2vector: Too many input arguments.');
    end

    try

        if ~iscell(value) && ~isstruct(value)
            vector = value;
            return
        end

        N = length(value);

        if nargin == 1
            % Return the values of a cell as a vector
            vector = cell2mat(value);
        else
            % Return the field of a struct from the cell as a vector
            try

                for i = N:-1:1
                    vector(:, i) = value{i}.(field);
                end

            catch
                vector = value.(field);
            end

        end

    catch
        vector = value{1};
    end

end
