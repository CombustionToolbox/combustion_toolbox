function self = complete_initialize(self, species)
    % Complete initialization process
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     species (cell): List of reactants
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases
    
    % If FLAG_INITIALIZE = false, then initialize self variable
    if self.Misc.FLAG_INITIALIZE
        return
    end

    % Get list of possible products
    if isempty(self.S.LS)
        LS = find_products(self, species);
    else
        LS = self.S.LS;
    end
    
    % Assign list of possible products
    self = list_species(self, LS);

    % Set Contained elements
    self = contained_elements(self);

    % Check species in DB, cataloged list of species, and compute stoichiometric matrix 
    self = initialize(self);
    
    % Update FLAG_INITIALIZE
    self.Misc.FLAG_INITIALIZE = true;
end