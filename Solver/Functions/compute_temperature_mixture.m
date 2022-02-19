function T = compute_temperature_mixture(self, species, moles, temperatures)
    % Compute equilibrium temperature [K] of a gaseous mixture compound of
    % n species with species at different temperatures.
    %
    % Inputs:
    %   - self: struct variable with all the data of the App (databases)
    %   - species: cell array with the species of the mixture
    %   - moles: vector or cell array with the number of moles of each species
    %   - temperatures: vector or cell array with the temperatures of each species
    % Output:
    %   - T: scalar with the temperature of the mixture at equilibrium
    
    try
        % Convert moles and temperature into vectors if needed
        moles = cell2vector(moles);
        temperatures = cell2vector(temperatures);
        if length(unique(temperatures)) > 1
            % Obtain specific heat capacity at constant volume and thermal enthalpy
            % by evaluating each species at its initial temperature
            cV_0 = get_property_list_species(@species_cV, species, temperatures, self.DB);
            DhT_0 = get_property_list_species(@species_DhT, species, temperatures, self.DB);
            % Guess of the temperature of the mixture at equilibrium
            T = sum(moles .* temperatures .* cV_0) / sum(moles .* cV_0); 
            % Constants
            it = 0; ERR = 1.0;
            % Solver: numerical 1D Newton-Raphson method
            while abs(ERR) > self.TN.tol0 && it < self.TN.itMax
                it = it + 1;
                % Obtain thermal enthalpy by evaluating each species at the guess
                % equilibrium temperature
                DhT = get_property_list_species(@species_DhT, species, T, self.DB);
                % Evaluate function
                FT = sum(moles .* DhT_0) - sum(moles .* DhT);
                % Perturb evaluated function
                DT = T*0.02;
                T_per = T + DT;
                DhT_per = get_property_list_species(@species_DhT, species, T_per, self.DB);
                FTX = sum(moles .* DhT) - sum(moles.*DhT_per);
                % Compute derivate (1D)
                DFTDT = (FTX - FT) / DT;
                % Solve equation to obtain the correction factor
                DeltaT = -DFTDT \ FT;
                % Apply correction
                T = T + DeltaT;
                % Compute error
                ERR = abs(DeltaT);
            end
        else
            % Same temperature
            T = temperatures(1);
        end
    catch
        T = max(temperatures);
    end
end

% SUB-PASS FUNCTIONS
function value = get_property_list_species(fun, species, temperatures, DB)
    % Get property evaluated at the given vector of temperatures and the
    % corresponded cell array of species
    try
        for i=length(species):-1:1
            value(i) = fun(species{i}, temperatures(i), DB); 
        end
    catch
        for i=length(species):-1:1
            value(i) = fun(species{i}, temperatures, DB); 
        end
    end
end