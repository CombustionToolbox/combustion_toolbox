function [FLAG_RE, REname] = isRefElm(reference_elements, species, T)
    % Check if the given species is a reference element
    %
    % Args:
    %     reference_elements (cell): List of reference elements with temperature intervals [K]
    %     species (str): Chemical species
    %     T (float): Temperature
    %
    % Return:
    %     name (str): Full name of the given species

    % Change lowercase 'l' to uppercase 'L' for Al, Cl, Tl, and Fl
    species(strfind(species, 'Al') + 1) = 'L';
    species(strfind(species, 'Cl') + 1) = 'L';
    species(strfind(species, 'Tl') + 1) = 'L';
    species(strfind(species, 'Fl') + 1) = 'L';

    FLAG_RE = false;
    REname = [];

    % Look for entries in the Reference_form_of_elements_with_T_intervals list
    % that partially match with the desired species and then check each one
    % sucessivelly

    j = find(contains(reference_elements, species));

    for i = 1:length(j)
        % disp(num2str(i))
        TentativeRefElm = reference_elements{j(i)};
        % Detect temperature interval
        n1 = strfind(TentativeRefElm, '[');
        n2 = strfind(TentativeRefElm, '-');
        n3 = strfind(TentativeRefElm, ']');
        T1 = sscanf(TentativeRefElm(n1 + 1:n2 - 1), '%f');
        T2 = sscanf(TentativeRefElm(n2 + 1:n3 - 1), '%f');

        if (T1 <= T) && (T <= T2)
            % Detect location of open parenthesis
            n_open_parenthesis = strfind(TentativeRefElm(1:n1 - 2), '(');
            % Detect location of '2'
            n_two = strfind(TentativeRefElm(1:n1 - 2), '2');
            % If thera are no '2's or parenthesis, the Species is essentially a noble gas
            if isempty(n_open_parenthesis) && isempty(n_two)
                % 1
                if strcmp(TentativeRefElm(1:n1 - 2), species)
                    FLAG_RE = true;
                    REname = TentativeRefElm(1:n1 - 2);
                end

            end

            % If there are '2's the species may be in the reference state or
            % not (e.g. O2 is, but O is not)
            if ~isempty(n_two)
                % 2
                if strcmp(TentativeRefElm(1:n_two), species)
                    FLAG_RE = true;
                    REname = TentativeRefElm(1:n_two);
                end

                if strcmp(TentativeRefElm(1:n_two - 1), species)
                    REname = TentativeRefElm(1:n_two);
                end

            end

            % If there are opening parenthesis, the species is in condensed phase
            if ~isempty(n_open_parenthesis)
                % 3
                if strcmp(TentativeRefElm(1:n_open_parenthesis - 1), species)
                    FLAG_RE = true;
                    REname = TentativeRefElm(1:n1 - 2);
                end

            end

        end

    end

end
