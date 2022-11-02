function append_cell = append_cells(cell1, cell2, varargin)
    % Append two or more cells in one common cell
    %
    % Args:
    %     cell1 (struct): Cell 1
    %     cell2 (struct): Cell 2
    %
    % Optional Args:
    %     celli (struct): Additional cells
    %
    % Returns:
    %     append_cell (struct): Merged cell

    % First case
    append_cell = fun_append(cell1, cell2);
    % Additional cases
    for i = 1:nargin - 2
        cell2 = varargin{i};
        append_cell = fun_append(append_cell, cell2);
    end

end

% SUB-PASS FUNCTIONS
function append_cell = fun_append(cell1, cell2)
    append_cell = [cell1, cell2];
end
