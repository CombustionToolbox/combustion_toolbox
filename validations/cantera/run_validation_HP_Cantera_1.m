function run_validation_HP_Cantera_1()
    % Run test validation_HP_Cantera_1:
    % Contrasted with: NASA's Chemical Equilibrium with Applications software
    % Problem type: Adiabatic T and composition at constant p
    % Temperature [K]   = 300;
    % Pressure    [bar] = 1;
    % Equivalence ratio [-] = 0.5:0.01:4
    % Initial mixture: C2H2_acetylene + AIR_IDEAL (79% N2 + 21% O2)
    % List of species considered: Cantera GRI-Mech 3.0

    % Define the gas mixture
    gas = Solution('gri30.yaml'); % Load GRI-Mech 3.0 mechanism
    graphite = Solution('graphite.yaml');

    % Initial conditions
    T0 = 300; % Initial temperature in K
    P0 = 1 * 1e5; % Initial pressure in Pa (1 bar)
    equivalence_ratios = 0.5:0.01:4; % Equivalence ratio range

    % Define the fuel and oxidizer
    fuel = 'C2H2'; % Acetylene
    oxidizer = {'O2', 'N2'};
    oxidizer_ratio = [1, 3.76]; % Air composition: 21% O2, 79% N2
    fuel_moles = 1;
    phi_st_fun = @(x, y, z) (x + y/4 - z/2);
    phi_st = phi_st_fun(2, 2, 0);

    % Initialize arrays to store results
    nSpecies = numel(speciesNames(gas)); % Get the number of species
    T_ad = zeros(size(equivalence_ratios));

    % Loop over each equivalence ratio
    for i = 1:length(equivalence_ratios)
        % Set the equivalence ratio
        phi = equivalence_ratios(i);
        
        % Oxidizer mole fractions (normalized by the total oxidizer ratio)
        oxidizer_moles = phi_st / phi .* oxidizer_ratio;
        
        % Total moles
        molesTotal = sum([fuel_moles, oxidizer_moles]);

        % Create a mole fraction array for the mixture
        X = zeros(1, nSpecies);
        X(speciesIndex(gas, fuel)) = fuel_moles / molesTotal;
        X(speciesIndex(gas, oxidizer{1})) = oxidizer_moles(1) / molesTotal;
        X(speciesIndex(gas, oxidizer{2})) = oxidizer_moles(2) / molesTotal;

        % Set the gas state
        set(gas, 'T', T0, 'P', P0, 'X', X);
        set(graphite, 'T', T0, 'P', P0);

        % Create a mixture of 1 mole of gas, and 0 moles of solid carbon
        mix = Mixture({gas, 1.0; graphite, 0});

        % Equilibrate the gas at constant pressure and enthalpy (HP problem)
        equilibrate(mix, 'HP');

        % Store the adiabatic flame temperature
        T_ad(i) = temperature(mix);
        
        % Store the mole fractions of all species
        % species_mole_fractions(i, :) = moleFractions(mix);
    end

    % Display results for some key species
    displaySpecies = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2',...
                      'HCN','H','OH','O','CN','NH3','CH4','C2H4','CH3',...
                      'NO','HCO','NH2','NH','N','CH'};

    % disp('Equivalence Ratio    Adiabatic Temperature (K)');
    % disp([equivalence_ratios' T_ad']);
    figure;
    plot(equivalence_ratios, T_ad);
    
    % for j = 1:length(displaySpecies)
    %     idx = speciesIndex(gas, displaySpecies{j});
    %     disp(['Species: ', displaySpecies{j}]);
    %     disp([equivalence_ratios' species_mole_fractions(:, idx)]);
    % end
end
