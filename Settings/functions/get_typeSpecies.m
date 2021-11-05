function typeSpecies = get_typeSpecies(app)
    % Create cell array with the type of species in the mixture
    typeFuel     = create_cell_ntimes('Fuel', length(app.PD.N_Fuel));
    typeOxidizer = create_cell_ntimes('Oxidizer', length(app.PD.N_Oxidizer));
    typeInert    = create_cell_ntimes('Inert', length(app.PD.N_Inert));
    typeSpecies  = [typeInert, typeOxidizer, typeFuel];
end