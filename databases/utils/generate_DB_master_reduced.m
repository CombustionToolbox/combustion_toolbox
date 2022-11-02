function DB_master_reduced = generate_DB_master_reduced(DB_master)
    % Generate Reduced Mater Database (DB_master_reduced) with the
    % thermodynamic data of the chemical species
    %
    % Args:
    %     DB_master (struct): Database with the thermodynamic data of the chemical species
    %
    % Returns:
    %     DB_master_reduced (struct): Reduced database with the thermodynamic data of the chemical species

    NameSpecies = fieldnames(DB_master);
    NSpecies = numel(NameSpecies);
    index = ones(NSpecies, 1);
    pattern = {'plus', 'minus', 'AL', 'Ag', 'F', 'CL', 'B', 'Ca', 'I', 'K', ...
                'Li', 'M', 'D', 'S', 'Rb', 'Pb', 'V', 'W', 'Z', 'G', 'T', 'Cd', 'Co', ...
                'Cr', 'Cs', 'Cu', 'Ni', 'U', 'Na', 'Nb', 'Hg', 'CP', 'HP'};

    for i = 1:NSpecies

        if contains(NameSpecies(i), pattern)
            index(i) = 0;
        elseif startsWith(NameSpecies(i), 'P')
            index(i) = 0;
        else
            DB_master_reduced.(NameSpecies{i}) = DB_master.(NameSpecies{i});
        end

    end

end
