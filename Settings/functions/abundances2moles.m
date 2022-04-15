function Xi = abundances2moles(elements, filename)
    % Read solar abundances in log 10 scale and compute the initial molar
    % fractions in the mixture [-]
    % 
    % Args:
    %
    %
    % Returns:
    %

    % Read abundances from filename
    [log_abundances, elements_file] = read_abundances(filename);
    % Reorganize abundances as the given element cell
    for i = length(elements):-1:1
        index(i) = find_ind(elements_file, elements{i});
    end
    % Compute molar fractions relative to H of the remaining elements in the mixture
    Xi = 10.^(log_abundances(index) - 12);
end