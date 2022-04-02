function Delta_n = compute_change_moles_gas_reaction(element_matrix, swtCondensed)
    % In order to compute the internal energy of formation from the enthalpy of
    % formation of a given species, we must determine the change in moles of
    % gases during the formation reaction of a mole of that species starting
    % from the elements in their reference state. 
    % 
    % Notes:
    % The only elements that are stable as diatomic gases are elements
    % 1 (H), 8 (N), 9 (O), 10 (F), and 18 (Cl). The remaining elements that
    % are stable as (monoatomic) gases are the noble gases He (3), Ne (11),
    % Ar (19), Kr (37), Xe (55), and Rn (87), which do not form any compound.
    %
    % Args: 
    %   element_matrix (float): Element matrix of the species
    %   swtCondensed (float):   0 or 1 indicating gas or condensed species
    %
    % Returns:
    %   Delta_n (float): Change in moles of gases during the formation reaction of a mole of that species starting from the elements in their reference state

    Delta_n_per_mole = sum(element_matrix(1,:)==[1, 8, 9, 10, 18]')/2 ... 
                     + sum(element_matrix(1,:)==[3, 11, 19, 37, 55, 87]');
    Delta_n = 1 - swtCondensed - dot(Delta_n_per_mole, element_matrix(2,:));
end