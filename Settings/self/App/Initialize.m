function self = Initialize(self)
    % Check if minors products species are contained in DB
    [self.DB, self.E, self.S, self.C] = check_DB(self, self.DB_master, self.DB);
    % Sort species: first gaseous species, secondly condensed species
    self = list_phase_species(self, self.S.LS);
    % Stoichiometric Matrix
    self = Stoich_Matrix(self);
end


%% SUB-PASS FUNCTIONS
function self = list_phase_species(self, LS)
    for ind=1:length(LS)
        Species = LS{ind};
        if ~self.DB.(Species).swtCondensed
           self.S.ind_nswt = [self.S.ind_nswt, ind];
        else
           self.S.ind_swt = [self.S.ind_swt, ind];
           if ~self.DB.(Species).ctTInt
              self.S.ind_cryogenic = [self.S.ind_cryogenic, ind];
           end
        end
    end
    self.S.ind_nswt = unique(self.S.ind_nswt);
    self.S.ind_swt  = unique(self.S.ind_swt);
    self.S.ind_cryogenic = unique(self.S.ind_cryogenic);
    self.S.LS = self.S.LS([self.S.ind_nswt, self.S.ind_swt]);
    self.S.NS = length(self.S.LS);
    self.S.NG = length(self.S.ind_nswt);
end


function self = Stoich_Matrix(self)
    self.C.A0.value = zeros(self.S.NS,self.E.NE);
    self.C.M0.value = zeros(self.S.NS,12);
    for i=1:self.S.NS
        txFormula = self.DB.(self.S.LS{i}).txFormula;
        self.DB.(self.S.LS{i}).Element_matrix = set_element_matrix(txFormula,self.E.elements);
        self.C.A0.value(i,self.DB.(self.S.LS{i}).Element_matrix(1,:)) = self.DB.(self.S.LS{i}).Element_matrix(2,:);
        self.C.M0.value(i,10) = self.DB.(self.S.LS{i}).swtCondensed;    
    end
    self.C.N0.value = self.C.M0.value(:, [1, 10]);
end