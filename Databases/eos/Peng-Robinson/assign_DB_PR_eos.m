function DB = assign_DB_PR_eos(DB)
    % Assing parameters for pure species considering Peng-Robinson Equation
    % of State (EoS)
    %
    %
    % References:
    % [1] THERMODYNAMIC AND THERMOPHYSICAL PROPERTIES OF HUMID AIR BY USING CUBIC PENG-ROBINSON EOS
    % [2] PROPERTIES OF GASES, Isidoro Martinez

    DB = add_species(DB, 'Ar', 151, 48.6 * 1e5, 0); % Ref. [2]

    DB = add_species(DB, 'CH4', 190.6, 46 * 1e5, 0.008);
    DB = add_species(DB, 'C2H4', 283, 51.2 * 1e5, 0.085); % Ref. [2]
    DB = add_species(DB, 'C2H6', 305.4, 48.84 * 1e5, 0.098);
    DB = add_species(DB, 'C3H8', 369.8, 42.46 * 1e5, 0.152);

    DB = add_species(DB, 'CO', 132.86, 34.94 * 1e5, 0.0497); % Ref. [1]
    DB = add_species(DB, 'CO2', 304.2, 73.76 * 1e5, 0.225);

    DB = add_species(DB, 'O2', 154.6, 50.4 * 1e5, 0.025); % Ref. [1]
    DB = add_species(DB, 'N2', 126.2, 33.9 * 1e5, 0.039); % Ref. [1]
    DB = add_species(DB, 'NO', 180, 64.8 * 1e5, 0.588); % Ref. [1]
    DB = add_species(DB, 'NO2', 261.9, 101 * 1e5, 0.834); % Ref. [1]

    DB = add_species(DB, 'H2', 33, 13.2 * 1e5, -0.22); % Ref. [2]
    DB = add_species(DB, 'H2O', 647.13, 22.34 * 1e5, 0.344);
    DB = add_species(DB, 'H2ObLb', 647.13, 22.34 * 1e5, 0.344);
end

% SUB-PASS FUNCTIONS
function DB = add_species(DB, species, T_c, p_c, w)
    % Coefficients
    DB.(species).critical_temperature = T_c; % [K]
    DB.(species).critical_pressure = p_c; % [Pa]
    DB.(species).acentric_factor = w; % [-]
end
