function max_rel_error_prop = run_test_Database_1()% value, listSpecies, DB
    % run_test_Database:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Equilibrium composition at defined T and p
    % Temperature [K]   = value;
    % Pressure    [bar] = 1.01325;
    % Initial mixture: Si + 9 C6H5OH_phenol
    % List of species considered: All (see routine find_products)
    
    % Import packages
    import combustiontoolbox.databases.NasaDatabase

    % Get NASA database
    DB = NasaDatabase();

    % 
    % % Get results Combustion Toolbox
    % 
    % % Load results Combustion Toolbox (validated)
    % 
    % % Load results CEA 
    % results_CEA = data_CEA(filename, display_species);
    % % Compute error
    % % * Molar fractions
    % max_rel_error_moles = compute_error_moles_CEA(results_CT, results_CEA, 'T', value, 'Xi', display_species);
    % % * Properties mixture 2
    % properties_y = {'cp', 'cv', 'h', 'e', 'g', 's'};
    % properties_x = create_cell_ntimes('T', length(properties_y));
    % max_rel_error_prop = compute_error_prop_CEA(results_CT, results_CEA, properties_x, value, properties_y, 'mix2');



    % Default values
    species = {'O2', 'N2', 'eminus', 'CO2'};
    properties = {'cp', 'h', 'e', 's', 'g'};
    T = 0:50:40000;

    % Unpack
    % for i = 1:2:nargin
    %     switch lower(varargin{i})
    %         case {'species', 'ls'}
    %             species = varargin{i + 1};
    %         case {'properties', 'prop', 'thermo'}
    %             properties = varargin{i + 1};
    %         case {'temperature', 't', 'temp'}
    %             T = varargin{i + 1};
    %     end
    % 
    % end
    
    % Definitions
    NS = length(species);
    NP = length(properties);
    NT = length(T);
    
    % Initialization
    Misc = Miscellaneous();
    
    % Compute thermodynamic polynomials
    [values_NASA, values_CT, y_labelname] = get_thermo(DB, species, properties, T, NS, NP, NT);
    
    % Get FLAG if the calculations imply extrapolation
    [FLAG_EXTRAPOLATION_PRE, FLAG_EXTRAPOLATION_POST] = get_FLAG_EXTRAPOLATION(DB, species, T, NS, NT);
    ALL_FLAGS = FLAG_EXTRAPOLATION_POST & FLAG_EXTRAPOLATION_PRE;
    
    % Define color palette
    color_palette = brewermap(NS, Misc.config.colorpalette);
    
    % Define normalized size of the figure
    Misc.config.innerposition = [0.1, 0.2, 0.8, 0.7];
    Misc.config.outerposition = [0.1, 0.2, 0.8, 0.7];
    
    % Plot results
    set_figure(Misc.config);
    
    tiledlayout(1, NP);
    
    for j = 1:NP
        ax = nexttile;
        Misc.config.xscale = 'log';
        Misc.config.yscale = 'log';
        ax = set_figure(ax, Misc.config);

        for i = 1:NS
            temp = reshape(values_CT(i, j, :), 1, NT);
            ax = plot_figure('T', T(~ALL_FLAGS(i, :)), y_labelname{j}, temp(~ALL_FLAGS(i, :)), 'color', color_palette(i, :), 'linestyle', '-', 'ax', ax);
            ax = plot_figure('T', T(ALL_FLAGS(i, :)), y_labelname{j}, temp(ALL_FLAGS(i, :)), 'color', color_palette(i, :), 'linestyle', ':', 'ax', ax);
            plot_figure('T', T, y_labelname{j}, reshape(values_NASA(i, j, :), 1, NT), 'color', color_palette(i, :), 'linestyle', '--', 'ax', ax);
        end
    
    end

end

% SUB-PASS FUNCTIONS
function [values_NASA, values_CT, y_labelname] = get_thermo(DB, species, properties, T, NS, NP, NT)
    % Compute thermodynamic polynomials
    for i = NS:-1:1
        for j = NP:-1:1
            for k = NT:-1:1
                % Get functions handles
                [funname_NASA, funname_CT, y_labelname{i}] = set_inputs_thermo_validations(properties{j});
    
                % Compute NASA
                values_NASA(i, j, k) = funname_NASA(species{i}, T(k), DB);
        
                % Compute CT
                values_CT(i, j, k) = funname_CT(species{i}, T(k), DB);
            end
    
        end
    
    end

end

function [FLAG_EXTRAPOLATION_PRE, FLAG_EXTRAPOLATION_POST] = get_FLAG_EXTRAPOLATION(DB, species, T, NS, NT)
    % Get FLAG if the calculations imply extrapolation 
    FLAG_EXTRAPOLATION_PRE = false(NS, NT);
    FLAG_EXTRAPOLATION_POST = false(NS, NT);
    for i = NS:-1:1
        FLAG_EXTRAPOLATION_PRE(i, :) = T < DB.(species{i}).T(1);
        FLAG_EXTRAPOLATION_POST(i, :) = T > DB.(species{i}).T(end);
    end

end

function [funname_NASA, funname_CT, y_labelname] = set_inputs_thermo_validations(property)
    % Set corresponding thermodynamic functions for NASA and Combustion
    % Toolbox
    %
    % Args:
    %     property (str): Thermodynamic property name
    %
    % Returns:
    %     Tuple containing
    %
    %     * funname_NASA (function): Function to use NASA's polynomials
    %     * funname_CT (function): Function to use Combustion Toolbox polynomials
    %     * y_labelname (str): Label y axis

    switch lower(property)
        case 'cp'
            funname_NASA = @species_cP_NASA;
            funname_CT = @species_cP;
            y_labelname = '\rm{Molar\ heat\ capacity\ at\ constant\ pressure\ [J/(mol-K)]}';
        case 'cv'
            funname_NASA = @species_cV_NASA;
            funname_CT = @species_cV;
            y_labelname = '\rm{Molar\ heat\ capacity\ at\ constant\ volume\ [J/(mol-K)]}';
        case {'det', 'et'}
            funname_NASA = @species_DeT_NASA;
            funname_CT = @species_DeT;
            y_labelname = '\rm{Molar\ thermal\ internal\ energy\ [kJ/mol]}';
        case {'dht', 'ht'}
            funname_NASA = @species_DhT_NASA;
            funname_CT = @species_DhT;
            y_labelname = '\rm{Molar\ thermal\ enthalpy\ [kJ/mol]}';
        case {'g0', 'g'}
            funname_NASA = @species_g0_NASA;
            funname_CT = @species_g0;
            y_labelname = '\rm{Molar\ Gibbs\ energy\ [kJ/mol]}';
        case {'h0', 'h'}
            funname_NASA = @species_h0_NASA;
            funname_CT = @species_h0;
            y_labelname = '\rm{Molar\ Enthalpy\ [kJ/mol]}';
        case {'s0', 's'}
            funname_NASA = @species_s0_NASA;
            funname_CT = @species_s0;
            y_labelname = '\rm{Molar\ Entropy\ [kJ/(mol-K)]}';
        case {'e', 'e0'}
            funname_NASA = @species_e0_NASA;
            funname_CT = @species_e0;
            y_labelname = '\rm{Internal\ energy\ [kJ/mol]}';
        case {'gamma', 'gammas'}
            funname_NASA = @species_gamma_NASA;
            funname_CT = @species_gamma;
            y_labelname = '\rm{Adiabatic\ index\ [-]}';
        otherwise
            error('There is not such property on files');
    end

end
