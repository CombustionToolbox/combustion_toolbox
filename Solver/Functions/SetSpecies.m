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

    for n = length(moles):-1:1
        mmi = self.DB.(species{n}).mm;                    % [g/mol]
        mi  = moles(n) * mmi * 1e-3;                      % [kg]
        hfi = self.DB.(species{n}).hf/1000;               % [kJ]
        efi = self.DB.(species{n}).ef/1000;               % [kJ]
        swtCondensed = self.DB.(species{n}).swtCondensed; % [bool]
        if length(self.DB.(species{n}).T) > 1
            DhTi = species_DhT(species{n},T,self.DB);     % [kJ]
            DeTi = species_DeT(species{n},T,self.DB);     % [kJ]
            cPi  = species_cP(species{n},T,self.DB);      % [J/K]
            cVi  = species_cV(species{n},T,self.DB);      % [J/K]
            s0i  = species_s0(species{n},T,self.DB);      % [kJ/K]
            if ~swtCondensed
                pVi = moles(n) * R0 * T * 1e-5; % For ideal gases [bar m3]
            else
                pVi = 0; % For condensed species [bar m3]
            end
        else
            DhTi = 0; % [kJ]
            DeTi = 0; % [kJ]
            cPi  = 0; % [J/K]
            cVi  = 0; % [J/K]
            s0i  = 0; % [kJ/K]
            pVi  = 0; % [bar m3]
        end   
        ind = find_ind(self.S.LS, species(n));
        M(ind, :) = [moles(n), moles(n) * [hfi, DhTi, efi, DeTi, cPi, cVi, s0i], pVi, swtCondensed, mi, mmi];
    end
end
