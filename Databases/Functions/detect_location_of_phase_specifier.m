function index_open_parenthesis = detect_location_of_phase_specifier(species) 
    % Detect the location of the opening pharenthesis of the phase identifier (if any) 
    %
    % Args:
    %     species (str): Chemical species
    %
    % Returns:
    %     n_open_parenthesis (float): Index of the location of the open parenthesis

    if ~isempty(find(species(:)=='(', 1))                                 % If the species name includes parenthesis
        index_open_parenthesis = find(species(:)=='(');                   % then detect the location of the open parenthesis
        index_close_parenthesis = find(species(:)==')');                  % and detect the location of the closing parenthesis
        if index_close_parenthesis(end) == length(species(:))             % If there is a closing parenthesis at the end of the species name, then it must be the closing parenthesis of the phase specifier
            index_open_parenthesis = index_open_parenthesis(end);         % then the phase specifier begins with the last parenthesys
        else                                                              % Some other times, the phase specifier is located before a short alphabetic name (e.g., 'C8H18(L),isooct') which is always preceeded by a comma 
            index_comma = find(species(:)==',');                          % detect the location of the comma (if any)
            if ~isempty(index_comma)                                      % If there is a comma (e.g., 'C8H18(L),isooct') then the last parenthesis may be the phase specifier as well, even though it is not located at the end of the species name
                if (index_comma-index_close_parenthesis(end)) == 1        % However, this is only true if the last parenthesis is just followed by the comma
                    index_open_parenthesis = index_open_parenthesis(end); % then the phase specifier begins also with the last parenthesys
                else
                    index_open_parenthesis = length(species)+1;           % Valid for species names w/o phase specifiers
                end
            else
                index_open_parenthesis = length(species)+1;               % Valid for species names w/o phase specifiers
            end
        end
    else
        index_open_parenthesis = length(species)+1;                       % Valid for species names w/o phase specifiers
    end
end