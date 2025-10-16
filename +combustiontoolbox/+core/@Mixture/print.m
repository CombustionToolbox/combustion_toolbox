function print(mix, varargin)
    % Print properties and composition of the given mixtures in the command
    % window
    %
    % Args:
    %     mix1 (Mixture): Mixture object with the properties of the mixture
    %
    % Optional Args:
    %     * mix2 (Mixture): Mixture object with the properties of the mixture
    %     * mixi (Mixture): Mixture object with the properties of the mixture
    %     * mixN (Mixture): Mixture object with the properties of the mixture
    %
    % Examples:
    %     * print(mix1)
    %     * print(mix1, mix2)
    %     * print(mix1, mix2, mix3)
    %     * print(mix1, mix2, mix3, mix4)
    
    % Definitions
    mintolDisplay = mix.config.mintolDisplay;
    compositionUnits = mix.config.compositionUnits;
    FLAG_COMPACT = mix.config.FLAG_COMPACT;

    % Unpack cell with mixtures
    mix = [{mix}, varargin(:)'];

    % Definitions
    numMixtures = nargin;
    listSpecies = mix{1}.chemicalSystem.listSpecies;
    problemType = mix{end}.problemType;
    
    % Check if problemType is specified
    if isempty(problemType), problemType = ''; end

    % Start
    fprintf('************************************************************************************************************\n');

    % Print header
    header_composition = print_header(problemType, numMixtures, mix);

    % Print properties
    print_properties(problemType, numMixtures, mix);
    
    % Print composition
    print_composition(mix, listSpecies, compositionUnits, header_composition, mintolDisplay);

    % End
    fprintf('************************************************************************************************************\n\n\n');

    % NESTED FUNCTIONS
    function print_composition(mix, listSpecies, units, header, mintolDisplay)
        % Print composition of the mixture in the command window
        %
        % Args:
        %     mix (struct): Struct with the properties of the mixture
        %     listSpecies (cell): Cell with the names of the species
        %     units (char): Units of the composition
        %     header (char): Header of the composition
        %     mintolDisplay (float): Minimum value to be displayed
    
        % Print composition (compact)
        if FLAG_COMPACT
            print_compact_composition(mix, listSpecies, units, mintolDisplay);
            return
        end
    
        % Print composition (sequential)
        for i = 1:numMixtures
            print_composition_sequential(mix{i}, listSpecies, compositionUnits, header{i}, mintolDisplay);
        end
    
    end

end

% SUB-PASS FUNCTIONS
function value = get_properties(property, numberMixtures, mix)
    % Get properties of the given mixtures
    %
    % Args:
    %     property (function/char): Function/char to get the property
    %     numberMixtures (float): Number of mixtures
    %     mix (cell): Cell with the properties of the N mixtures
    %
    % Returns:
    %     value (float): Value of the property
    %
    % Examples:
    %     * get_properties(@temperature, 1, mix)
    %     * get_properties('uShock', 1, mix)

    if ischar(property)
        value = cell2vector(mix, property);
        return
    end

    for i = numberMixtures:-1:1
        value(i) = property(mix{i});
    end

end

function line = set_string_value(Nmixtures, varargin)
    % Set the char to print the properties
    %
    % Args:
    %     Nmixtures (float): Number of mixtures
    %
    % Optional Args:
    %     * format (char): Format to print the properties
    %     * limiter (char): Limiter to print the properties
    %
    % Returns:
    %     line (char): Char to print the properties
    
    % Default
    format = '%12.4f';
    limiter = '|';

    % Unpack additional inputs
    if nargin > 1
        format = varargin{1};
        limiter = varargin{2};
    end

    % Set line
    line_body = sprintf('   %s  %s', format, limiter);
    line_end = sprintf('   %s\n', format);
    line = [repmat(line_body, 1, Nmixtures - 1), line_end];
end

function print_properties(ProblemType, numberMixtures, mix)
    % Print properties of the mixture in the command window
    %
    % Args:
    %     ProblemType (char): Type of problem
    %     numberMixtures (float): Number of mixtures
    %     mix (cell): Cell with the properties of the N mixtures

    % Definitions
    string_value = set_string_value(numberMixtures);
    string_value_2 = set_string_value(numberMixtures - 1);
    
    % Print properties
    fprintf(['T [K]          |', string_value], get_properties(@temperature, numberMixtures, mix));
    fprintf(['p [bar]        |', string_value], get_properties(@pressure, numberMixtures, mix));
    fprintf(['r [kg/m3]      |', string_value], get_properties(@density, numberMixtures, mix));
    fprintf(['h [kJ/kg]      |', string_value], get_properties(@enthalpy_mass, numberMixtures, mix));
    fprintf(['e [kJ/kg]      |', string_value], get_properties(@intEnergy_mass, numberMixtures, mix));
    fprintf(['g [kJ/kg]      |', string_value], get_properties(@gibbs_mass, numberMixtures, mix));
    fprintf(['s [kJ/(kg-K)]  |', string_value], get_properties(@entropy_mass, numberMixtures, mix));
    fprintf(['W [g/mol]      |', string_value], get_properties(@MolecularWeight, numberMixtures, mix));
    fprintf(['(dlV/dlp)T [-] |', string_value], get_properties('dVdp_T', numberMixtures, mix));
    fprintf(['(dlV/dlT)p [-] |', string_value], get_properties('dVdT_p', numberMixtures, mix));
    fprintf(['cp [kJ/(kg-K)] |', string_value], get_properties(@cp_mass, numberMixtures, mix));
    fprintf(['gamma [-]      |', string_value], get_properties(@adiabaticIndex, numberMixtures, mix));
    fprintf(['gamma_s [-]    |', string_value], get_properties(@adiabaticIndex_sound, numberMixtures, mix));
    fprintf(['sound vel [m/s]|', string_value], get_properties(@soundspeed, numberMixtures, mix));

    if contains(ProblemType, 'SHOCK') || contains(ProblemType, 'DET')
        fprintf(['u [m/s]        |', string_value], [get_properties(@velocity_relative, 1, mix(1)), get_properties('uShock', numberMixtures - 1, mix(2:end))]);
        fprintf(['Mach number [-]|', string_value], [get_properties(@velocity_relative, 1, mix(1)), get_properties('uShock', numberMixtures - 1, mix(2:end))] ./ get_properties(@soundspeed, numberMixtures, mix));
    end

    if contains(ProblemType, '_OBLIQUE') || contains(ProblemType, '_POLAR')
        fprintf('------------------------------------------------------------------------------------------------------------\n');
        fprintf('PARAMETERS\n');
        fprintf(['min wave  [deg]|                 |', string_value_2], get_properties('betaMin', numberMixtures - 1, mix(2:end)));
        fprintf(['wave angle[deg]|                 |', string_value_2], get_properties('beta', numberMixtures - 1, mix(2:end)));
        fprintf(['deflection[deg]|                 |', string_value_2], get_properties('theta', numberMixtures - 1, mix(2:end)));

        if contains(ProblemType, '_POLAR')
            fprintf(['max def.  [deg]|                 |', string_value_2], get_properties('thetaMax', numberMixtures - 1, mix(2:end)));
            fprintf(['sonic def.[deg]|                 |', string_value_2], get_properties('thetaSonic', numberMixtures - 1, mix(2:end)));
        end

    elseif contains(ProblemType, 'PRANDTL_MEYER')
        fprintf('------------------------------------------------------------------------------------------------------------\n');
        fprintf('PARAMETERS\n');
        fprintf(['deflection[deg]|                 |', string_value_2], get_properties('theta', numberMixtures - 1, mix(2:end)));
    elseif contains(ProblemType, 'ROCKET')
        string_value_3 = set_string_value(numberMixtures - 2);

        fprintf('------------------------------------------------------------------------------------------------------------\n');
        fprintf('PERFORMANCE PARAMETERS\n');
        fprintf(['A/At [-]       |                 |                 |', string_value_3], get_properties('areaRatio', numberMixtures - 2, mix(3:end)));
        fprintf(['CSTAR [m/s]    |                 |                 |', string_value_3], get_properties('cstar', numberMixtures - 2, mix(3:end)));
        fprintf(['CF [-]         |                 |                 |', string_value_3], get_properties('cf', numberMixtures - 2, mix(3:end)));
        fprintf(['Ivac [s]       |                 |                 |', string_value_3], get_properties('I_vac', numberMixtures - 2, mix(3:end)));
        fprintf(['Isp  [s]       |                 |                 |', string_value_3], get_properties('I_sp', numberMixtures - 2, mix(3:end)));
    end

    fprintf('------------------------------------------------------------------------------------------------------------\n');
end

function print_compact_composition(mixCell, listSpecies, units, mintolDisplay)
    % mixCell         = cell array of mixture objects
    % listSpecies     = cell array of species names
    % units           = 'molar fraction', 'mass fraction', or 'mol'
    % mintolDisplay  = threshold below which species are considered minor

    % Definitions
    numMixtures = numel(mixCell);
    nSpecies    = numel(listSpecies);

    % Build the composition matrix: each column is one mixture's composition
    comp_matrix = zeros(nSpecies, numMixtures);
    for m = 1:numMixtures
        switch lower(units)
            case 'mol'
                comp_matrix(:, m) = moles(mixCell{m});
                short_label = 'Ni [mol]';
            case 'molar fraction'
                comp_matrix(:, m) = moleFractions(mixCell{m});
                short_label = 'Xi [-]';
            case 'mass fraction'
                comp_matrix(:, m) = massFractions(mixCell{m});
                short_label = 'Yi [-]';
            otherwise
                error('Unsupported composition unit: %s', units);
        end
    end

    % Determine which species exceed mintolDisplay in ANY mixture (major species)
    major_mask  = any(comp_matrix > mintolDisplay, 2);
    major_vals  = comp_matrix(major_mask, :);
    major_names = listSpecies(major_mask);

    % Sort major species by their composition in the FIRST mixture (descending)
    [~, idxSort] = sort(major_vals(:,1), 'descend');
    major_vals   = major_vals(idxSort, :);
    major_names  = major_names(idxSort);

    % Print composition property
    fprintf('COMPOSITION    ', short_label);
    for m = 1:numMixtures
        fprintf('%15s   ', short_label);
    end
    fprintf('\n');

    % Prepare a single format string for aligned column
    line = set_string_value(numMixtures, '%12.4e', ' ');

    % Print each major species in one row, columns for each mixture
    for i = 1:size(major_vals, 1)
        fprintf('%-16s', major_names{i});  % Species name, left-justified
        fprintf(line, major_vals(i,:));    % Species composition in each mixture
    end

    % Compute MINORS (sum of all species below threshold in EVERY mixture)
    minor_mask = ~major_mask; % Get the mask for minor species
    Nminor = sum(minor_mask); % Number of minor species
    minor_values = sum(comp_matrix(minor_mask, :), 1); % Sum of minor species for each mixture

    % Print a single MINORS row with columns for each mixture
    fprintf('%-16s', sprintf('MINORS[+%d]', Nminor));
    fprintf([line, '\n'], minor_values);  

    % Print TOTAL row for each mixture
    totals = sum(comp_matrix, 1);  % sum of all species for each mixture
    fprintf('%-16s', 'TOTAL');
    fprintf(line, totals);
end


function print_composition_sequential(mix, listSpecies, units, header, mintolDisplay)
    % Print composition of the mixture in the command window
    %
    % Args:
    %     mix (struct): Struct with the properties of the mixture
    %     listSpecies (cell): Cell with the names of the species
    %     units (char): Units of the composition
    %     header (char): Header of the composition
    %     mintolDisplay (float): Minimum value to be displayed
    
    switch lower(units)
        case 'mol'
            variable = moles(mix);
            short_label = 'Ni [mol]\n';
        case 'molar fraction'
            variable = moleFractions(mix);
            short_label = '  Xi [-]\n';
        case 'mass fraction'
            variable = massFractions(mix);
            short_label = '  Yi [-]\n';
    end

    fprintf([header, short_label]);
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [variable, ind_sort] = sort(variable, 'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = variable > mintolDisplay;
    minor = sum(variable(~j));

    for i = 1:length(j)

        if j(i)
            fprintf('%-20s %1.4e\n', listSpecies{ind_sort(i)}, variable(i));
        end

    end

    Nminor = length(variable) - sum(j);
    s_space_Nminor = char(32 * ones(1, 4 - numel(num2str(Nminor))));
    fprintf('MINORS[+%d] %s     %12.4e\n\n', Nminor, s_space_Nminor, minor);
    fprintf('TOTAL            %14.4e\n', sum(variable));
    fprintf('------------------------------------------------------------------------------------------------------------\n');
end

function header_composition = print_header(problemType, numberMixtures, mix)
    % Print header in the command window
    %
    % Args:
    %     problemType (char): Type of the problem
    %     numberMixtures (float): Number of mixtures
    %     mix (cell): Cell with the properties of the N mixtures
    %
    % Returns:
    %     header_composition (cell): Cell with the header indicating the type/state of the mixture
    
    FLAG_PHI = ~isempty(mix{1}.equivalenceRatio);

    if problemType & FLAG_PHI
        fprintf('------------------------------------------------------------------------------------------------------------\n');
        fprintf('Problem type: %s  | Equivalence ratio = %4.4f\n', problemType, mix{1}.equivalenceRatio);
    elseif problemType & ~FLAG_PHI
        fprintf('------------------------------------------------------------------------------------------------------------\n');
        fprintf('Problem type: %s\n', problemType);
    elseif FLAG_PHI
        fprintf('------------------------------------------------------------------------------------------------------------\n');
        fprintf('Equivalence ratio = %4.4f\n', mix{1}.equivalenceRatio);
    end

    fprintf('------------------------------------------------------------------------------------------------------------\n');

    if contains(problemType, '_OBLIQUE') || contains(problemType, '_POLAR')

        if numberMixtures == 2
            header_composition = {'STATE 1               ', 'STATE 2               '};
            fprintf('               |     STATE 1     |     STATE 2\n');
        elseif numberMixtures == 3
            header_composition = {'STATE 1               ', 'STATE 2-W             ', 'STATE 2-S             '};
            fprintf('               |     STATE 1     |     STATE 2-W   |     STATE 2-S\n');
        elseif numberMixtures == 4
            header_composition = {'STATE 1               ', 'STATE 2               ', 'STATE 3-W             ', 'STATE 3-S             '};
            fprintf('               |     STATE 1     |     STATE 2     |     STATE 3-W   |     STATE 3-S\n');
        end

    elseif contains(problemType, '_R')
        header_body = '     STATE %d     |';
        header_end = '     STATE %d\n';
        header = [repmat(header_body, 1, numberMixtures - 1), header_end];
        value = 1:numberMixtures;
        fprintf(['               |', header], value);

        header_composition = {'STATE 1               ', ...
                              'STATE 2               ', ...
                              'STATE 3               ', ...
                              'STATE 4               ', ...
                              'STATE 5               ', ...
                              'STATE 6               ', ...
                              'STATE 7               ', ...
                              'STATE 8               '};
    elseif contains(problemType, 'ROCKET')

        if numberMixtures == 3
            header_composition = {'INLET CHAMBER         ', ...
                                  'OUTLET CHAMBER        ', ...
                                  'THROAT                '};
            fprintf('               |  INLET CHAMBER  | OUTLET CHAMBER  |      THROAT \n');
        elseif numberMixtures > 3

            if numberMixtures > 4
                header_exit_prop_last = '|      EXIT\n';
            else
                header_exit_prop_last = [];
            end

            if mix{3}.areaRatio == 1
                header_exit = repmat({'EXIT                  '}, 1, numberMixtures - 3);
                header_composition = {'INLET CHAMBER         ', ...
                                      'OUTLET CHAMBER        ', ...
                                      'THROAT                ', ...
                                    header_exit{1:end}};

                if numberMixtures > 4
                    header_exit_prop = [repmat({'|      EXIT       |'}, 1, numberMixtures - 5), header_exit_prop_last];
                    fprintf(['               |  INLET CHAMBER  | OUTLET CHAMBER  |     THROAT      ', header_exit_prop{1:end}]);
                else
                    fprintf('               |  INLET CHAMBER  | OUTLET CHAMBER  |     THROAT      |      EXIT\n');
                end

            else
                header_exit = repmat({'EXIT                  '}, 1, numberMixtures - 4);
                header_composition = {'INLET CHAMBER         ', ...
                                      'INJECTOR              ', ...
                                      'OUTLET CHAMBER        ', ...
                                      'THROAT                ', ...
                                    header_exit{1:end}};

                if numberMixtures > 4
                    header_exit_prop = [repmat({'|      EXIT       |'}, 1, numberMixtures - 5), header_exit_prop_last];
                    fprintf(['               |  INLET CHAMBER  |     INJECTOR    | OUTLET CHAMBER  |     THROAT      ', header_exit_prop{1:end}]);
                else
                    fprintf('               |  INLET CHAMBER  |     INJECTOR    | OUTLET CHAMBER  |     THROAT \n');
                end

            end

        end

    elseif problemType & numberMixtures == 2
        header_composition = {'REACTANTS             ', ...
                        'PRODUCTS              '};
        fprintf('               |    REACTANTS    |      PRODUCTS\n');
    else
        header_body = '    MIXTURE %d    |';
        header_end = '    MIXTURE %d\n';
        header = [repmat(header_body, 1, numberMixtures - 1), header_end];
        value = 1:numberMixtures;
        fprintf(['               |', header], value);

        header_composition = {'MIXTURE 1             ', ...
                              'MIXTURE 2             ', ...
                              'MIXTURE 3             ', ...
                              'MIXTURE 4             ', ...
                              'MIXTURE 5             ', ...
                              'MIXTURE 6             ', ...
                              'MIXTURE 7             ', ...
                              'MIXTURE 8             '};
    end

end
