function gui_ReactantsValueChanged(self, event)
    % Update reactants & GUI
    self = gui_update_Reactants(self, event);
    % Update equivalence ratio 
    self.edit_phi2.Value = self.edit_phi.Value;
end

% SUB-PASS FUNCTIONS
function self = gui_update_Reactants(self, event)
    FLAG_IDEAL_AIR = true;
    % Get temperature of the mixture
    T = sscanf(self.PR1.Value, '%f');
    switch self.Reactants.Value
        case '1' % No species selected
            gui_empty_Reactants(self);
            return
        case '2' % AIR
            FLAG_FUEL = false;
            phi_ratio = 1;
            [~, ~, Data] = compute_reactants(phi_ratio, T, FLAG_IDEAL_AIR, FLAG_FUEL);
            self.UITable_R.Data = Data;
        case '3' % METHANE + AIR
            FLAG_FUEL = true;
            phi_ratio = 2;
            [Ni, Xi, Data] = compute_reactants(phi_ratio, T, FLAG_IDEAL_AIR, FLAG_FUEL);
            self.UITable_R.Data = [Data; {'CH4', Ni(3), Xi(3), 'Fuel', T}];
        otherwise % SET NEW SPECIES
            try
                species = seeker_species(self, event);
            catch
                message = {'Species not found.'};
                uialert(self.UIFigure, message, 'Warning', 'Icon', 'warning');
                return
            end
            % Add species to the table
            self.UITable_R.Data = [self.UITable_R.Data;{species, 0, 0, 'Fuel', T}];
    end
    self.UITable_P.Data = self.UITable_R.Data(:, 1);    % Update UITable_P  (species, numer of moles, mole fractions, temperature)
    self.UITable_R2.Data = self.UITable_R.Data(:, 1:3); % Update UITable_R2 (species, numer of moles, mole fractions)
    % Find if there are fuels in the mixture to compute the equivalence ratio
    self = gui_get_typeSpecies(self);
    % Compute equivalence ratio
    
end

function self = gui_empty_Reactants(self)
    self.UITable_R.Data  = [];
    self.UITable_P.Data  = [];
    self.UITable_R2.Data = [];
end

function species = seeker_species(self, event)
    seekSpecies = event.Value;
    listSpecies = self.S.LS_DB;
    index = 0;
    seekIndex = false;
    while ~seekIndex
        index = index + 1;
        seekIndex = startsWith(listSpecies{index}, seekSpecies, 'IgnoreCase', false);
    end
    if index
        species = listSpecies{index}; % Species found in the Database
    else
        species = []; % Species not found in the Database
    end
end

function [Ni, Xi, Data] = compute_reactants(phi_ratio, Temperature, FLAG_IDEAL_AIR, FLAG_FUEL)
    Ni = compute_N_air(phi_ratio, FLAG_IDEAL_AIR);
    if FLAG_FUEL
        Ni = [Ni, 1];
    end
    Xi = Ni / sum(Ni);
    if FLAG_IDEAL_AIR
        Data = {'N2', Ni(1), Xi(1), 'Inert', Temperature; 'O2', Ni(2), Xi(2), 'Oxidant', Temperature};
    else
        Data = {'N2', Ni(1), Xi(1), 'Inert', Temperature; 'O2', Ni(2), Xi(2), 'Oxidant', Temperature; 'Ar', Ni(3), Xi(3), 'Inert', Temperature; 'CO2', Ni(4), Xi(4), 'Oxidant', Temperature};
    end
end

function Ni_air = compute_N_air(phi_ratio, FLAG_IDEAL_AIR)
    if FLAG_IDEAL_AIR
        Ni_air = [phi_ratio * 0.79/0.21, phi_ratio];
    else
        Ni_air = [phi_ratio * 0.78084/0.209476, phi_ratio, phi_ratio * 0.009365/0.209476, phi_ratio * 0.009365/0.209476];
    end
end