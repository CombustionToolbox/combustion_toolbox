function self = define_FOI(self, i)
    % Set up mixture: fuel, oxidizer and diluent/inert species
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     i (float): Position of the evaluated problem
    %
    % Returns:
    %     self (struct): Data of the mixture, conditions, and databases

    species = [self.PD.S_Fuel, self.PD.S_Oxidizer, self.PD.S_Inert];
    % Check reactant species are contained in the list of products (initial computations)
    self = check_FOI(self, species);
    % Define moles Fuel
    if ~self.Misc.FLAG_N_Fuel && ~self.Misc.FLAG_GUI
        self.PD.N_Fuel = 1;
    end

    % Computation of theoretical equivalence ratio
    self = define_F(self);
    % Define moles Oxidizer
    if ~self.Misc.FLAG_N_Oxidizer && ~isempty(self.PD.S_Oxidizer)

        if self.Misc.FLAG_WEIGHT
            self.PD.N_Oxidizer = self.PD.phi_t / self.PD.phi.value(i) .* convert_weight_percentage_to_moles(self.PD.S_Oxidizer, self.PD.wt_ratio_oxidizers, self.DB) / sum(self.PD.N_Fuel);
            %             self.PD.N_Oxidizer = sum(self.PD.N_Fuel) ./ (self.PD.phi.value(i) .* self.PD.phi_t .* sum(convert_weight_percentage_to_moles(self.PD.S_Oxidizer, self.PD.wt_ratio_oxidizers, self.DB)));
        else
            self.PD.N_Oxidizer = self.PD.phi_t / self.PD.phi.value(i) .* self.PD.ratio_oxidizers_O2;
        end

    end

    % Define moles Diluent/Inert
    if ~self.Misc.FLAG_N_Inert && ~isempty(self.PD.S_Inert)

        if self.Misc.FLAG_WEIGHT
            self.PD.N_Inert = self.PD.phi_t / self.PD.phi.value(i) .* convert_weight_percentage_to_moles(self.PD.S_Inert, self.PD.wt_ratio_inerts, self.DB);
        else
            self.PD.N_Inert = self.PD.phi_t / self.PD.phi.value(i) .* self.PD.ratio_inerts_O2;
        end

    end

    % Set inputs
    moles = [self.PD.N_Fuel, self.PD.N_Oxidizer, self.PD.N_Inert];
    temperatures = [self.PD.T_Fuel, self.PD.T_Oxidizer, self.PD.T_Inert];
    % Compute Temperature of the reactants (species with different temperature)
    self.PD.TR.value = compute_temperature_mixture(self, species, moles, temperatures);
    % Define mixture at equilibrium temperature
    self = define_F(self);
    self = define_O(self);
    self = define_I(self);
    % Compute property matrix of the reactives for the given conditions
    R = self.PD.R_Fuel + self.PD.R_Oxidizer + self.PD.R_Inert;
    R(:, 1:self.C.M0.ind_ni - 1) = self.C.M0.value(:, 1:self.C.M0.ind_ni - 1);
    % Compute properties of the reactives for a given temperature and
    % pressure
    self.PS.strR{i} = compute_properties(self, R, self.PD.pR.value, self.PD.TR.value);
    self.PS.strR{i}.phi = self.PD.phi.value(i);
    self.PS.strR{i}.LS = merged_cells({self.PD.S_Fuel, self.PD.S_Oxidizer, self.PD.S_Inert});
    [~, ind_LS, ~] = intersect(self.S.LS, self.PS.strR{i}.LS);
    [~, ind_LS_Inert, ~] = intersect(self.S.LS, self.PD.S_Inert);
    self.PS.strR{i}.ind_LS = ind_LS;
    self.PS.strR{i}.ind_LS_Inert = ind_LS_Inert;
    self.PS.strR{i}.phi_c = compute_phi_c(self.PD.Fuel);
    % Compute percentage Fuel, Oxidizer/Fuel ratio and equivalence ratio
    self = compute_ratios_fuel_oxidizer(self, R, i);
end

% SUB-PASS FUNCTIONS
function merged = merged_cells(cells)
    merged = [];

    for i = 1:length(cells)
        merged = [merged, cells{i}];
    end

end

function mass = compute_mass_mixture(self, properties_matrix)
    % Compute mass mixture [kg]

    % Compute total number of moles [mol]
    N = sum(properties_matrix(:, self.C.M0.ind_ni));
    % Compute mean molecular weight [g/mol]
    W = dot(properties_matrix(:, self.C.M0.ind_ni), properties_matrix(:, self.C.M0.ind_W)) / N;
    % Compute mass mixture [kg]
    mass = W * N * 1e-3; % [kg]
end

function self = compute_ratios_fuel_oxidizer(self, R, i)
    % Compute percentage Fuel, Oxidizer/Fuel ratio and equivalence ratio
    if ~isempty(self.PD.S_Fuel) && ~isempty(self.PD.S_Oxidizer)
        mass_fuel = compute_mass_mixture(self, self.PD.R_Fuel);
        mass_oxidizer = compute_mass_mixture(self, self.PD.R_Oxidizer);
        mass_mixture = compute_mass_mixture(self, R);
        self.PS.strR{i}.percentage_Fuel = mass_fuel / mass_mixture * 100;
        self.PS.strR{i}.OF = mass_oxidizer / mass_fuel;
        self.PS.strR{i}.FO = 1 / self.PS.strR{i}.OF;
        self.PS.strR{i}.FO_moles = sum(self.PD.R_Fuel(:, self.C.M0.ind_ni)) / sum(self.PD.R_Oxidizer(self.S.ind_ox_ref, self.C.M0.ind_ni));
        self.PS.strR{i}.FO_moles_st = abs(sum(self.PD.R_Fuel(:, self.C.M0.ind_ni)) / (self.PS.strR_Fuel.x + self.PS.strR_Fuel.x2 + self.PS.strR_Fuel.x3 + self.PS.strR_Fuel.y / 4 - self.PS.strR_Fuel.z / 2) * (0.5 * self.S.atoms_ox_ref));
        self.PS.strR{i}.phi = self.PS.strR{i}.FO_moles / self.PS.strR{i}.FO_moles_st;
    elseif ~isempty(self.PD.S_Fuel)
        self.PS.strR{i}.percentage_Fuel = 100;
        self.PS.strR{i}.FO = inf;
        self.PS.strR{i}.OF = 0;
        self.PS.strR{i}.phi = '-';
    else
        self.PS.strR{i}.percentage_Fuel = 0;
        self.PS.strR{i}.FO = 0;
        self.PS.strR{i}.OF = inf;
        self.PS.strR{i}.phi = '-';
    end

end