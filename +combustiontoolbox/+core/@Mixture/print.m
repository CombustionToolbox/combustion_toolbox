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
    %     * printMixture(mix1)
    %     * printMixture(mix1, mix2)
    %     * printMixture(mix1, mix2, mix3)
    %     * printMixture(mix1, mix2, mix3, mix4)
    
    % Temporal
    mintol_display = 1e-14; % (will be moved to Miscellaneous)
    composition_units = 'molar fraction'; % Possible values: mol, molar fraction or mass fraction
    
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
    for i = 1:numMixtures
        print_composition(mix{i}, listSpecies, composition_units, header_composition{i}, mintol_display);
    end

    % End
    fprintf('************************************************************************************************************\n\n\n');
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

function line = set_string_value(Nmixtures)
    % Set the char to print the properties
    %
    % Args:
    %     Nmixtures (float): Number of mixtures
    %
    % Returns:
    %     line (char): Char to print the properties
    
    line_body = '   %12.4f  |';
    line_end = '   %12.4f\n';
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

function print_composition(mix, LS, units, header, mintol_display)
    % Print composition of the mixture in the command window
    %
    % Args:
    %     mix (struct): Struct with the properties of the mixture
    %     LS (cell): Cell with the names of the species
    %     units (char): Units of the composition
    %     header (char): Header of the composition
    %     mintol_display (float): Minimum value to be displayed
    
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
    j = variable > mintol_display;
    minor = sum(variable(~j));

    for i = 1:length(j)

        if j(i)
            fprintf('%-20s %1.4e\n', LS{ind_sort(i)}, variable(i));
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
