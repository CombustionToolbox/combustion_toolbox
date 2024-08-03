function titlename = getTitle(obj)
    % Get a title based on the problem type and species involved
    %
    % Args:
    %     obj (Mixture): Mixture object
    %
    % Returns:
    %     titlename (char): Title based on the problem type and species involved

    % Definitions
    FLAG_FUEL = ~isempty(obj.listSpeciesFuel);
    FLAG_OXIDIZER = ~isempty(obj.listSpeciesOxidizer);
    FLAG_INERT = ~isempty(obj.listSpeciesInert);
    FLAG_RATIO_INERTS_O2 = ~isempty(obj.ratioOxidizer);
    N_oxidizer = length(obj.listSpeciesOxidizer);
    molesFuel = obj.molesFuel;
    listSpeciesFuel = obj.listSpeciesFuel;
    
    % Initialization
    label_problemtype = [];

    % Check if the type of species are specified
    if ~FLAG_FUEL && ~FLAG_OXIDIZER && ~FLAG_INERT
        FLAG_PASS = obj.quantity > 0;
        molesFuel = obj.quantity(FLAG_PASS);
        listSpeciesFuel = obj.listSpecies(FLAG_PASS);
    end

    % Get the title based on the problem type and species involved
    if ~isempty(obj.problemType)
        label_problemtype = [strrep(obj.problemType, '_', ' '), ': '];
    end

    titlename = [label_problemtype, cat_moles_species(molesFuel, listSpeciesFuel)];

    if ~FLAG_OXIDIZER && ~FLAG_INERT
        return
    end

    if FLAG_FUEL
        titlename = [titlename, ' + $\frac{', sprintf('%.3g', obj.stoichiometricMoles), '}{\phi}$'];
    end

    if FLAG_OXIDIZER

        if N_oxidizer > 1 && FLAG_FUEL
            titlename = [titlename, '('];
        end

        ind = combustiontoolbox.utils.findIndex(obj.listSpeciesOxidizer, 'O2');

        if ind
            obj.molesOxidizer = obj.molesOxidizer / obj.molesOxidizer(ind);
        end

        titlename = [titlename, cat_moles_species(obj.ratioOxidizer, obj.listSpeciesOxidizer)];
    end

    if FLAG_INERT && FLAG_RATIO_INERTS_O2
        titlename = [titlename, ' + ', cat_moles_species(obj.molesInert, obj.listSpeciesInert)];
    end

    if FLAG_OXIDIZER && N_oxidizer > 1 && FLAG_FUEL
        titlename = [titlename, ')'];
    end

    if FLAG_INERT && ~FLAG_RATIO_INERTS_O2
        titlename = [titlename, ' + ', cat_moles_species(obj.molesInert, obj.listSpeciesInert)];
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
    
    cat_text = strrep(cat_text, '_$_{', '$_{');
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

    cat_text = [value, combustiontoolbox.utils.display.species2latex(species)];
end