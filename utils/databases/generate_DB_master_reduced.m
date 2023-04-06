function DB_master_reduced = generate_DB_master_reduced(DB_master)
    % Generate Reduced Mater Database (DB_master_reduced) with the
    % thermodynamic data of the chemical species (deprecated)
    %
    % Args:
    %     DB_master (struct): Database with the thermodynamic data of the chemical species
    %
    % Returns:
    %     DB_master_reduced (struct): Reduced database with the thermodynamic data of the chemical species

    species = fieldnames(DB_master);
    NS = numel(species);
    index = ones(NS, 1);
    pattern = {'plus', 'minus', 'AL', 'Ag', 'F', 'CL', 'B', 'Ca', 'I', 'K', ...
               'Li', 'M', 'D', 'S', 'Rb', 'Pb', 'V', 'W', 'Z', 'G', 'T', 'Cd', 'Co', ...
               'Cr', 'Cs', 'Cu', 'Ni', 'U', 'Na', 'Nb', 'Hg', 'CP', 'HP'};

    for i = 1:NS

        if contains(species(i), pattern)
            index(i) = 0;
        elseif startsWith(species(i), 'P')
            index(i) = 0;
        else
            DB_master_reduced.(species{i}) = DB_master.(species{i});
        end

    end

end
