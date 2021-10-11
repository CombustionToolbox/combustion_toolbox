function n_open_parenthesis = detect_location_of_phase_specifier(Species) 

% Detect the location of the opening pharentesis of the phase identifier (if any) 

if ~isempty(find(Species(:)=='(', 1))                              % If the species name includes parenthesis
    n_open_parenthesis = find(Species(:)=='(');                 % then detect the location of the open parenthesis
    n_close_parenthesis = find(Species(:)==')');                %  and detect the location of the closing parenthesis
    if n_close_parenthesis(end) == length(Species(:))           % If there is a closing parenthesis at the end of the species name, then it must be the closing parenthesis of the phase specifier
        n_open_parenthesis = n_open_parenthesis(end);           % then the phase specifier begins with the last parenthesys
    else                                                        % Some other times, the phase specifier is located before a short alphabetic name (e.g., 'C8H18(L),isooct') which is always preceeded by a comma 
        n_comma = find(Species(:)==',');                        % Detect the location of the comma (if any)
        if ~isempty(n_comma)                                    % If there is a comma (e.g., 'C8H18(L),isooct') then the last parenthesis may be the phase specifier as well, even though it is not located at the end of the species name
            if (n_comma-n_close_parenthesis(end)) == 1          % however, this is only true if the last parenthesis is just followed by the comma
                n_open_parenthesis = n_open_parenthesis(end);   % then the phase specifier begins also with the last parenthesys
            else
                n_open_parenthesis = length(Species)+1;         % valid for species names w/o phase specifiers
            end
        else
            n_open_parenthesis = length(Species)+1;             % valid for species names w/o phase specifiers
        end
    end
else
    n_open_parenthesis = length(Species)+1;                     % valid for species names w/o phase specifiers
end