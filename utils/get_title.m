function titlename = get_title(self)
    % Get a title based on the problem type and species involved
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %
    % Returns:
    %     titlename (char): Title based on the problem type and species involved

    label_problemtype = strrep(self.PD.ProblemType, '_', ' ');
    titlename = [label_problemtype, ': ', cat_moles_species(self.PD.N_Fuel, self.PD.S_Fuel)];

    if isempty(self.PD.S_Oxidizer) && isempty(self.PD.S_Inert)
        return
    end

    if ~isempty(self.PD.S_Fuel)
        titlename = [titlename, ' + '];
    end

    titlename = [titlename, '$\ \frac{', sprintf('%.3g', self.PD.phi_t), '}{\phi}$'];

    if ~isempty(self.PD.S_Oxidizer)

        if length(self.PD.S_Oxidizer) > 1
            titlename = [titlename, '('];
        end

        ind = find_ind(self.PD.S_Oxidizer, 'O2');

        if ind
            self.PD.N_Oxidizer = self.PD.N_Oxidizer / self.PD.N_Oxidizer(ind);
        end

        titlename = [titlename, cat_moles_species(self.PD.N_Oxidizer, self.PD.S_Oxidizer)];
    end

    if ~isempty(self.PD.S_Inert) && ~isempty(self.PD.ratio_inerts_O2)
        titlename = [titlename, ' + ', cat_moles_species(self.PD.N_Inert, self.PD.S_Inert)];
    end

    if ~isempty(self.PD.S_Oxidizer) && length(self.PD.S_Oxidizer) > 1
        titlename = [titlename, ')'];
    end

    if ~isempty(self.PD.S_Inert) && isempty(self.PD.ratio_inerts_O2)
        titlename = [titlename, ' + ', cat_moles_species(self.PD.N_Inert, self.PD.S_Inert)];
    end

end

% SUB-PASS FUNCTIONS
function cat_text = cat_moles_species(moles, species)
    % This function takes an array of moles and a cell array of species
    % names and returns a string with the chemical formula for the mixture
    %
    % Args:
    %     moles (float): Vector containing the number of moles of each species
    %     species (cell): Cell containing the name of each species
    %
    % Returns:
    %     cat_text (char): Chemical formula for the mixture

    N = length(species);
    cat_text = [];

    if ~N
        return
    end

    cat_text = cat_mol_species(moles(1), species{1});

    for i = 2:N
        cat_text = [cat_text, ' + ', cat_mol_species(moles(i), species{i})];
    end

end

function cat_text = cat_mol_species(mol, species)
    % This function takes a number of moles and a species name and returns
    % a string with the chemical formula for the species
    %
    % Args:
    %     moles (float): Vector containing the number of moles of each species
    %     species (cell): Cell containing the name of each species
    %
    % Returns:
    %     cat_text (char): Chemical formula for the mixture

    if mol == 1
        value = [];
    else
        value = sprintf('%.3g', mol);
    end

    cat_text = [value, species2latex(species)];
end