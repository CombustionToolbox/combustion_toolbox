function cell = assign_vector2cell(cell, vector, varargin)
    % Assign values of a vector into a cell
    %
    % Args:
    %     cell (cell): Cell in which the values of the given vector are going to be included
    %     vector (any): Vector with the values that are going to be included in the cell
    %
    % Optional Args:
    %     ind (float): List of index positions to assign specific positions to the cell
    %
    % Returns:
    %     cell (cell): Cell with the values of the given vector

    if nargin > 2
        ind = varargin{1};

        if islogical(ind)
            % Convert to index positions
            ind = find(ind);
        end

        for i = length(vector):-1:1
            cell(ind(i)) = {vector(i)};
        end

    else

        for i = length(vector):-1:1
            cell(i) = {vector(i)};
        end

    end

end
