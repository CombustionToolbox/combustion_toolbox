function element_matrix = set_element_matrix(txFormula, elements)
    % Compute element matrix of the given species
    %
    % Args:
    %     txFormula (str): Chemical formula
    %
    % Returns:
    %     element_matrix(float): Element matrix
    %
    % Example:
    %     For CO2
    %
    %     element_matrix = [7, 9; 1, 2]
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

        element_i = txFormula(start0 - 2:end0 - 6);

        if strcmp(element_i, '  ')
        else

            if strcmp(element_i(2), ' ')
                element_i = element_i(1);
            end

            element_matrix(1, i) = find(strcmpi(elements, element_i));
            element_matrix(2, i) = sscanf(txFormula(start0:end0), '%f');
        end

    end

end
