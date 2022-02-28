function vector = cell2vector(varargin)
    % Convert values of a individual cell into a vector. If the value 
    % correspond with a struct it can return as a vector the values of a 
    % given fieldname.
    if nargin > 0
        str = varargin{1};
    end
    if nargin > 1
        field = varargin{2};
    end
    if nargin > 2
        error('Error in @cell2vector: Too many input arguments.');
    end
    
    try
        Nstruct = length(str);
        if nargin == 1
            % Return the values of a cell as a vector
            for i=Nstruct:-1:1
                vector(:, i) = str{i};
            end
        else
            % Return the field of a struct as a vector
            try
                for i=Nstruct:-1:1
                    vector(:, i) = str{i}.(field);
                end
            catch
                vector = str.(field);
            end
        end
    catch
        vector = varargin{1};
    end
end