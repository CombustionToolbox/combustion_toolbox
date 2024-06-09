function gui_update_phi(app)
    % Update GUI: equivalence ratio, O/F, and percentage Fuel
    %
    % Args:
    %     app (object): Combustion Toolbox app object
    %     self (object): Data of the mixture, conditions, and databases
    
    if isempty(app.mixture.equivalenceRatio)
        return
    end
    
    app.edit_phi.Value = sprintf('%.5g', round(app.mixture.equivalenceRatio, 5));
    app.edit_phi2.Value = app.edit_phi.Value;
    app.edit_phi3.Value = app.edit_phi.Value;
    app.edit_OF.Value = app.mixture.oxidizerFuelMassRatio;
    app.edit_F.Value = app.mixture.percentageFuel;
end