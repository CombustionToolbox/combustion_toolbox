function moles = abundances2moles(elements, filename, varargin)
    % Read solar abundances in log 10 scale and compute the initial molar
    % fractions in the mixture [-]
    %
    % Args:
    %    elements (cell): List with the given elements
    %    filename (file): Filename with the data
    %
    % Optional Args:
    %    metallicity (float): Metallicity
    %
    % Returns:
    %   moles (float): moles relative to H of the remaining elements in the mixture

    % Unpack
    metallicity = unpack(varargin);
    % Read abundances from filename assuming unity metallicity
    [log_abundances, elements_file] = read_abundances(filename);
    % Recompute with metallicity. NOTE: H and He do not change
    if nargin > 2
        index_H = find_ind(elements_file, 'H');
        index_He = find_ind(elements_file, 'He');
        index_change = 1:1:length(elements_file);
        index_change([index_H, index_He]) = [];
        log_abundances(index_change) = log_abundances(index_change) + log10(metallicity);
    end

    % Reorganize abundances as the given element cell
    for i = length(elements):-1:1
        index(i) = find_ind(elements_file, elements{i});
    end

    % Compute moles relative to H of the remaining elements in the mixture
    moles = 10.^(log_abundances(index) - 12);
end

% SUB-PASS FUNCTIONS
function metallicity = unpack(variable)
    % Unpack extra inputs
    metallicity = 1;

    if ~isempty(variable)
        metallicity = variable{1};
    end

end
