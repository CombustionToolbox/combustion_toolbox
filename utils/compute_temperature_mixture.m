function T = compute_temperature_mixture(problemType, database, listSpecies, moles, temperatures)
    % Compute equilibrium temperature [K] of a gaseous mixture compound of
    % n species with species at different temperatures
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     species (cell): Cell array with the species of the mixture
    %     moles (float): Moles of the species in the mixture [mol]
    %     temperatures (float): Vector or cell array with the temperatures of each species
    %
    % Returns:
    %     T (float): Temperature of the mixture at equilibrium
    %
    % Example:
    %     T = compute_temperature_mixture(self, {'H2', 'O2'}, [1, 2], [300, 400])
    
    % Definitions
    tol0 = 1e-3;
    itMax = 100;

    try
        % Convert moles and temperature into vectors if needed
        moles = cell2vector(moles(:)');
        temperatures = cell2vector(temperatures(:)');

        if isscalar(unique(temperatures))
            % Same temperature
            T = temperatures(1);
            return
        end

        switch lower(problemType)
            case {'tv', 'ev', 'sv', 'volume', 'v'}
                funCpOrCv = @species_cV;
                funHorE = @species_e0;
            otherwise
                funCpOrCv = @species_cP;
                funHorE = @species_h0;
        end

        % Obtain specific heat capacity at constant volume and thermal enthalpy
        % by evaluating each species at its initial temperature
        CpOrCv_0 = getPropertyListSpecies(funCpOrCv, listSpecies, temperatures, database);
        HorE_0 = getPropertyListSpecies(funHorE, listSpecies, temperatures, database);

        % Guess of the temperature of the mixture at equilibrium
        T = sum(moles .* temperatures .* CpOrCv_0) / sum(moles .* CpOrCv_0);

        % Initialization
        it = 0; STOP = 1.0;

        % Solver: numerical 1D Newton-Raphson method
        while STOP > tol0 && it < itMax
            % Update iterations
            it = it + 1;

            % Obtain enthalpy [J/mol] by evaluating each species at the guess
            % equilibrium temperature
            HorE = getPropertyListSpecies(funHorE, listSpecies, T, database);

            % Evaluate function
            f = sum(moles .* HorE_0) - sum(moles .* HorE);
            f_rel = f / sum(moles .* HorE);

            % Compute first derivative [J/(mol-K)]
            CpOrCv = getPropertyListSpecies(funCpOrCv, listSpecies, T, database);
            df = -sum(moles .* CpOrCv);

            % Compute correction factor
            DeltaT = -f / df;

            % Apply correction
            T = T + DeltaT;

            % Compute error
            STOP = max(abs([DeltaT, f_rel]));
        end

    catch

        if isempty(temperatures)
            T = obj.PD.TR.value;
        else
            T = max(temperatures);
        end

    end

end

% SUB-PASS FUNCTIONS
function value = getPropertyListSpecies(fun, listSpecies, temperatures, database)
    % Get property evaluated at the given vector of temperatures and the
    % corresponded cell array of species
    
    for i = length(listSpecies):-1:1
        try
            value(i) = fun(listSpecies{i}, temperatures(i), database);
        catch
            value(i) = fun(listSpecies{i}, temperatures, database);
        end
    end

end
