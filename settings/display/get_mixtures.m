function mixtures = get_mixtures(PS, pattern)
    % Get all non-empty mixture
    %
    % Args:
    %       PS (struct): Struct with all the data of Problem Solution (PS)
    %       pattern (str): Pattern/s name of the mixture
    %
    % Returns:
    %       mixtures (cell): Cell with all the non-empty mixtures
    %
    % Example:
    %
    %       mixtures = get_mixtures(self.PS, 'mix');
    %       mixtures = get_mixtures(self.PS, 'strP');

    % Get fieldnames
    fieldnames_value = fieldnames(PS);
    % Get all non-empty mixture for the given pattern
    mixtures = get_mix(PS, pattern, fieldnames_value);
end

% SUB-PASS FUNCTIONS
function mixtures = get_mix(PS, pattern, fieldnames_value)
    % Get all non-empty mixture for the given pattern

    % Initialization
    j = 0;
    mixtures = {};
    fieldnames_mix = fieldnames_value(contains(fieldnames_value, pattern));
    % Get non-empty mixtures
    for i = length(fieldnames_mix):-1:1

        if ~isempty(PS.(fieldnames_mix{i}))
            j = j + 1;
            mixtures{j} = PS.(fieldnames_mix{i});
        end

    end

end
