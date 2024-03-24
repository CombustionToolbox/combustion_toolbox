function elementMatrix = set_element_matrix(formula, elements)
    % Compute element matrix of the given species
    %
    % Args:
    %     formula (char): Chemical formula
    %
    % Returns:
    %     elementMatrix(float): Element matrix
    %
    % Example:
    %     For CO2
    %
    %     elementMatrix = [7, 9; 1, 2]
    %
    %     That is, the species contains 1 atom of element 7 (C) and
    %     2 atoms of element 9 (O)

    % Definitions
    N = 40;
    NE = 5;
    step = 8;
    
    % Fill element matrix
    for i = 5:-1:1
        end0 = N - step * (NE - i);
        start0 = end0 - 5;

        element_i = formula(start0 - 2:end0 - 6);

        if strcmp(element_i, '  ')
        else

            if strcmp(element_i(2), ' ')
                element_i = element_i(1);
            end

            elementMatrix(1, i) = find(strcmpi(elements, element_i));
            elementMatrix(2, i) = sscanf(formula(start0:end0), '%f');
        end

    end

end
