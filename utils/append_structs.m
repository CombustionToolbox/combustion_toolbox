function append_s = append_structs(s1, s2, varargin)
    % Append two or more structs in one common struct
    %
    % Args:
    %     s1 (struct): Struct 1
    %     s2 (struct): Struct 2
    %
    % Optional Args:
    %     si (struct): Additional structs
    %
    % Returns:
    %     append_s (struct): Merged struct

    % First case
    append_s = fun_append(s1, s2);
    % Additional cases
    for i = 1:nargin - 2
        s2 = varargin{i};
        append_s = fun_append(append_s, s2);
    end

end

% SUB-PASS FUNCTIONS
function append_s = fun_append(s1, s2)
    append_s = struct();
    fields = fieldnames(s1);
    N = length(fields);
    cat_cell = [struct2cell(s1), struct2cell(s2)];

    for i = 1:N
        append_s.(fields{i}) = cat(2, cat_cell{i, :});
    end

end
