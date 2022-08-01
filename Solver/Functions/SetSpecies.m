function M = SetSpecies(self, species, moles, T)
    % Fill the properties matrix with the data of the mixture
    %
    % Args:
    %     self (struct):  Data of the mixture, conditions, and databases
    %     species (cell): Species contained in the system
    %     moles (float):  Moles of the species in the mixture [mol]
    %     T (float):      Temperature [K]
    %
    % Returns:
    %     M (float):      properties matrix

    M = self.C.M0.value;
    R0 = self.C.R0;

    for i = length(moles):-1:1
        W = self.DB.(species{i}).mm; % [g/mol]
        mi  = moles(i) * W * 1e-3;   % [kg]
        hfi = self.DB.(species{i}).hf/1000; % [kJ/mol]
        efi = self.DB.(species{i}).ef/1000; % [kJ/mol]
        phase = self.DB.(species{i}).phase; % [bool]
        if length(self.DB.(species{i}).T) > 1
            h0i = species_h0(species{i},T,self.DB); % [kJ/mol]
            cPi = species_cP(species{i},T,self.DB); % [J/mol-K]
            s0i = species_s0(species{i},T,self.DB); % [kJ/mol-K]
            if ~phase
                pVi = moles(i) * R0 * T * 1e-5; % For ideal gases [bar m3]
            else
                pVi = 0; % For condensed species [bar m3]
            end
        else
            h0i = hfi; % [kJ/mol]
            cPi = 0; % [J/mol-K]
            s0i = 0; % [kJ/mol-K]
            pVi = 0; % [bar m3]
        end   
        ind = find_ind(self.S.LS, species(i));
        M(ind, :) = [moles(i), moles(i) * [hfi, h0i, efi, cPi, s0i], pVi, phase, mi, W];
    end
end
