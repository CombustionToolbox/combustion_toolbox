function struct_var = append_fields_struct(struct_var)

    fields = fieldnames(struct_var);
    N = length(fields);
    temp_cases = 0;
    temp_LS = [];

    for i = 1:N
        N_cases = numel(struct_var.(fields{i}));
        struct_var.(fields{i}) = reshape(struct_var.(fields{i}), N_cases, 1);
        temp_cases = max(temp_cases, N_cases);

        if strcmpi(fields{i}, 'X')

            for j = 1:length(struct_var.X)
                temp_LS = [temp_LS; struct_var.X(j).mole(:, 1)];
            end

            temp_LS = unique(temp_LS, 'stable');
            temp_X = zeros(length(temp_LS), temp_cases);
            % Fill matrix
            t = 1;

            for j = 1:length(struct_var.X)
                k = length(struct_var.X(j).mole{1, 2});
                index_current = find_ind(struct_var.X(j).mole(:, 1), temp_LS);
                index = find_ind(temp_LS, struct_var.X(j).mole(:, 1));

                for w = 1:length(index)
                    value = struct_var.X(j).mole{index_current(w), 2};

                    if isempty(value)
                        temp_X(index(w), t:(k + t - 1)) = 0;
                    else
                        temp_X(index(w), t:(k + t - 1)) = value;
                    end

                end

                t = k + 1;
            end

        end

    end

    struct_var.X = temp_X;
    struct_var.LS = temp_LS;
end
