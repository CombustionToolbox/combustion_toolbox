function self = initialize(self)
    % This routine has three tasks:
    %   - Check that all species are contained in the Database
    %   - Establish cataloged list of species according to the state of the
    %     phase (gaseous or condensed). It also obtains the indices of
    %     cryogenic liquid species, e.g., liquified gases
    %   - Compute Stoichiometric Matrix
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases

    try
        % Check if minors products species are contained in DB
        [self.DB, self.E, self.S, self.C] = check_DB(self, self.DB_master, self.DB);
        % Sort species: first gaseous species, secondly condensed species
        self = list_phase_species(self, self.S.LS);
        % Set Stoichiometric and Properties Matrix
        self = stoich_prop_matrix(self);
    catch ME
        errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
            ME.stack(1).name, ME.stack(1).line, ME.message);
        errorSolution = sprintf('Error Solution:\n%s', 'None.');
        fprintf('%s\n%s\n', errorMessage, errorSolution);
    end

end
