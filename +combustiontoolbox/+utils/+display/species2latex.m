function speciesLatex = species2latex(species, varargin)
    % Convert species name into LateX format
    %
    % Args:
    %     species (char): Species name
    %
    % Optional Args:
    %     FLAG_BURCAT (bool): If true, do not remove Burcat database prefix (default: true)
    %
    % Returns:
    %     speciesLatex (char): Species name in LateX format
    %
    % Examples:
    %     * species2latex('H2ObLb')
    %     * species2latex('Si2H6_M')
    %     * species2latex('Si2H6_M', false)
    
    % Default
    FLAG_BURCAT = true;

    % Unpack
    if nargin > 1
        FLAG_BURCAT = varargin{1};
    end

    % Remove database Millennium prefix '_M'
    if ~FLAG_BURCAT
        species = strrep(species, '_M', '');
    else
        species = strrep(species, '_M', '$_{\rm M}$');
    end

    % Check numbers
    index = regexp(species, '[0-9]');
    N = length(index);

    if ~N
        speciesLatex = species;
    else
        speciesLatex = [];
        pos1 = 1;

        for i = 1:N
            pos2 = index(i) - 1;
            speciesLatex = strcat(speciesLatex, species(pos1:pos2), '$_{', species(pos2 + 1), '}$');
            pos1 = pos2 + 2;
        end

        if pos1 < length(species) + 1

            if isletter(species(pos1)) || species(pos1) == '$'
                aux = 0;
            else
                aux = 1;
            end

            speciesLatex = [speciesLatex, species(pos1 + aux:end)];
        end

    end

    % Check if ions
    speciesLatex = strrep(speciesLatex, 'plus', '$^+$');
    speciesLatex = strrep(speciesLatex, 'minus', '$^-$');

    % Check if is there are parenthesis
    [ind_start, ind_end] = regexp(speciesLatex, '(?=b)(.*?)(?=b)');
    if ind_start
        speciesLatex(ind_start(1)) = '(';
        speciesLatex(ind_end(end) + 1) = ')';
    end
    
    % Check concatenate $$
    speciesLatex = strrep(speciesLatex, '$$', '');

    % Check suffix
    index = regexp(speciesLatex, '_[a-zA-Z]');
    speciesLatex(index) = '';
    
    % Joining the Burcat subscript with the previous subscript
    [~, ind_end] = regexp(speciesLatex, '(?=_{)(.*?)(?=})');
    FLAG_DOUBLE_SUBSCRIPTS = contains(speciesLatex, '_{\rm M}') & length(ind_end) > 1;

    if FLAG_DOUBLE_SUBSCRIPTS
        % Remove _M subscript
        speciesLatex = strrep(speciesLatex, '_{\rm M}', '');
        % Join _M subscript with the first one
        speciesLatex = [speciesLatex(1:ind_end(end-1)), ', \rm{M}', speciesLatex(ind_end(end-1) + 1:end)];
        % Check concatenate $$
        speciesLatex = strrep(speciesLatex, '$$', '');
    end

    % Check long numbers
    speciesLatex = strrep(speciesLatex, '}_{', '');
end