function excel_cell = get_excel_cell(mix, species, phi)
    % Construct a cell with the thermodynamic data of the given mixture
    %
    % Args:
    %     mix (struct): Properties of the mixture
    %     species (cell): List of species
    %     phi (float): Vector of equivalence ratio
    %
    % Returns:
    %     excell_cell (cell): Cell with the thermodynamic data of the given mixture
    % 
    % Example:
    %     excell_cell(self.PS.strR, self.S.LS, self.PD.phi.value)
    
    % Definitions
    N = length(phi);
    NP = 19;
    NS = numel(species);
    composition_unit = 'kg';
    composition_variable = 'X';
    FLAG_VELOCITY = isfield(mix{1}, 'u');
    FLAG_THERMO_DERIVATIVES = isfield(mix{1}, 'dVdT_p');

    % Get data phases species
    phase = mix{1}.phase';
    
    % Set phase labels
    phase_labels = cell(1, NS);
    phase_labels(phase > 0) = {'condensed'};
    phase_labels(phase == 0) = {'gas'};
    
    % Headers - thermodynamic properties
    excel_cell{1, 1} = 'phi';      excel_cell{2, 1} = '[-]';
    excel_cell{1, 2} = 'T';        excel_cell{2, 2} = '[K]';
    excel_cell{1, 3} = 'p';        excel_cell{2, 3} = '[bar]';
    excel_cell{1, 4} = 'rho';      excel_cell{2, 4} = '[kg/m^3]';
    excel_cell{1, 5} = 'h';        excel_cell{2, 5} = ['[kJ/', composition_unit,']'];
    excel_cell{1, 6} = 'e';        excel_cell{2, 6} = ['[kJ/', composition_unit,']'];
    excel_cell{1, 7} = 's';        excel_cell{2, 7} = ['[kJ/(', composition_unit,'·K)]'];
    excel_cell{1, 8} = 'g';        excel_cell{2, 8} = ['[kJ/', composition_unit,']'];
    excel_cell{1, 9} = 'N';        excel_cell{2, 9} = '[mol]';
    excel_cell{1, 10} = 'W';       excel_cell{2, 10} = '[g/mol]';
    excel_cell{1, 11} = 'MW';      excel_cell{2, 11} = '[g/mol]';
    excel_cell{1, 12} = 'm';       excel_cell{2, 12} = '[kg]';
    excel_cell{1, 13} = 'dVdT_p';  excel_cell{2, 13} = '[-]';
    excel_cell{1, 14} = 'dVdp_T';  excel_cell{2, 14} = '[-]';
    excel_cell{1, 15} = 'cP';      excel_cell{2, 15} = ['[kJ/(', composition_unit,'·K)]'];
    excel_cell{1, 16} = 'gamma';   excel_cell{2, 16} = '[-]';
    excel_cell{1, 17} = 'gamma_s'; excel_cell{2, 17} = '[-]';
    excel_cell{1, 18} = 'sound';   excel_cell{2, 18} = '[m/s]';
    excel_cell{1, 19} = 'u';       excel_cell{2, 19} = '[m/s]';
    
    % Headers - chemical species
    for j = NS:-1:1
        excel_cell{1, NP + j} = sprintf('%s_%s (%s)', composition_variable, species{j}, phase_labels{j});
        excel_cell{2, NP + j} = '[-]';
    end
    
    % Assign data
    for i = N:-1:1
        % Equivalence ratio
        excel_cell{i + 2, 1} = phi(i);
        % Temperature
        excel_cell{i + 2, 2} = temperature(mix{i});
        % Pressure
        excel_cell{i + 2, 3} = pressure(mix{i});
        % Density
        excel_cell{i + 2, 4} = density(mix{i});
        % Specific enthalpy
        excel_cell{i + 2, 5} = enthalpy_mass(mix{i});
        % Specific internal energy
        excel_cell{i + 2, 6} = intEnergy_mass(mix{i});
        % Specific Gibbs energy
        excel_cell{i + 2, 7} = gibbs_mass(mix{i});
        % Specific entropy
        excel_cell{i + 2, 8} = entropy_mass(mix{i});
        % Total number of moles
        excel_cell{i + 2, 9} = mix{i}.N;
        % Molecular weight (gases)
        excel_cell{i + 2, 10} = MolecularWeight(mix{i});
        % Mean molecular weight
        excel_cell{i + 2, 11} = meanMolecularWeight(mix{i});
        % Mass
        excel_cell{i + 2, 12} = mass(mix{i});

        % Thermodynamic derivatives
        if FLAG_THERMO_DERIVATIVES
            % Thermodynamic derivative at constant pressure
            excel_cell{i + 2, 13} = mix{i}.dVdT_p;
            % Thermodynamic derivative at constant temperature
            excel_cell{i + 2, 14} = mix{i}.dVdp_T;
        end

        % Specific heat at constant pressure
        excel_cell{i + 2, 15} = cp_mass(mix{i});
        % Specific heat ratio
        excel_cell{i + 2, 16} = adiabaticIndex(mix{i});
        % Adiabatic index
        excel_cell{i + 2, 17} = adiabaticIndex_sound(mix{i});
        % Sound velocity
        excel_cell{i + 2, 18} = soundspeed(mix{i});

        % Shock velocity
        if FLAG_VELOCITY
            excel_cell{i + 2, 19} = velocity_relative(mix{i});
        end
    
        % Molar fractions
        for j = NS:-1:1
            excel_cell{i + 2, NP + j} = mix{i}.Xi(j);
        end
        
    end

end