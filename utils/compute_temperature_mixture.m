function T = compute_temperature_mixture(self, species, moles, temperatures)
    % Compute equilibrium temperature [K] of a gaseous mixture compound of
    % n species with species at different temperatures.
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     species (str): Cell array with the species of the mixture
    %     moles: Moles of the species in the mixture [mol]
    %     temperatures: Vector or cell array with the temperatures of each species
    %
    % Returns:
    %     T (float): Temperature of the mixture at equilibrium

    try
        % Convert moles and temperature into vectors if needed
        moles = cell2vector(moles);
        temperatures = cell2vector(temperatures);

        if length(unique(temperatures)) > 1
            % Obtain specific heat capacity at constant volume and thermal enthalpy
            % by evaluating each species at its initial temperature
            cV_0 = get_property_list_species(@species_cV_NASA, species, temperatures, self.DB);
            h0_0 = get_property_list_species(@species_h0_NASA, species, temperatures, self.DB);
            % Guess of the temperature of the mixture at equilibrium
            T = sum(moles .* temperatures .* cV_0) / sum(moles .* cV_0);
            % Constants
            it = 0; ERR = 1.0;
            % Solver: numerical 1D Newton-Raphson method
            while abs(ERR) > self.TN.tol0 && it < self.TN.itMax
                it = it + 1;
                % Obtain enthalpy [kJ/mol] by evaluating each species at the guess
                % equilibrium temperature
                h0 = get_property_list_species(@species_h0_NASA, species, T, self.DB);
                % Evaluate function
                f = sum(moles .* h0_0) - sum(moles .* h0);
                f_rel = f / sum(moles .* h0);
                % Compute first derivative [kJ/(mol-K)]
                cP = get_property_list_species(@species_cP_NASA, species, T, self.DB) * 1e-3;
                df =- sum(moles .* cP);
                % Compute correction factor
                DeltaT =- f / df;
                % Apply correction
                T = T + DeltaT;
                % Compute error
                ERR = max(abs(DeltaT), f_rel);
            end

        else
            % Same temperature
            T = temperatures(1);
        end

    catch

        if isempty(temperatures)
            T = self.PD.TR.value;
        else
            T = max(temperatures);
        end

    end

end

% SUB-PASS FUNCTIONS
function value = get_property_list_species(fun, species, temperatures, DB)
    % Get property evaluated at the given vector of temperatures and the
    % corresponded cell array of species
    try

        for i = length(species):-1:1
            value(i) = fun(species{i}, temperatures(i), DB);
        end

    catch

        for i = length(species):-1:1
            value(i) = fun(species{i}, temperatures, DB);
        end

    end

end
