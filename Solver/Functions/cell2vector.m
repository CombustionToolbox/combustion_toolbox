function vector = cell2vector(value, varargin)
    % Convert values of a individual cell into a vector. If the value 
    % correspond with a struct it can return as a vector the values of a 
    % given fieldname.
    %
    % Args:
    %     value (cell/struct): Data of the mixture, conditions, and databases
    % Optional Args:
    %     field (str):         Fieldname of the given value (struct)
    %
    % Returns:
    %     vector (any):        Vector with the values of the individual cell/fieldname (struct)
    
    if nargin > 1
        field = varargin{1};
    end
    if nargin > 2
        error('Error in @cell2vector: Too many input arguments.');
    end
    
    try
        Nstruct = length(value);
        if nargin == 1
            % Return the values of a cell as a vector
            for i=Nstruct:-1:1
                vector(:, i) = value{i};
            end
        else
            % Return the field of a struct as a vector
            try
                for i=Nstruct:-1:1
                    vector(:, i) = value{i}.(field);
                end
            catch
                vector = value.(field);
            end
        end
    catch
        vector = value;
    end
end