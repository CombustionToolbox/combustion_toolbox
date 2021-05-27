function self = Initialize(self)
    % List of species that we are going to compute
    self = Compute_Species(self);
    % Index gaseous and condensed species
    self = list_phase_species(self, self.S.LS);
    % Sort species: first gaseous species, secondly condensed species
    self = rearrange_species(self);
    % Index fixed species
    self = Index_fixed_Species(self);
    % Stoichiometric Matrix
    self = Stoich_Matrix(self);
    % Compute CHON equilibria for the minors products considered
    self = Compute_minors_species(self);
    % CH4 major species
    self = CH4_major(self);
    % OH major species
    self = OH_major(self);
    % Ask problem type
    self.PD.ProblemType = Ask_problem(self);
end
%% NESTED FUNCTIONS
function self = Compute_Species(self)
    % First we eliminate from the minor species list those considered major
    % species in case the user has included any
    n_pass = [];
    for n = 1:self.M.Lminors
        if ~any(strcmp(self.M.minors_products{n}, self.S.LS_fixed))
            n_pass = [n_pass, n];
        end
    end
    self.M.minors_products = self.M.minors_products(n_pass);
    self.S.LS = [self.S.LS_fixed, self.M.minors_products];
    self.S.NS = length(self.S.LS);
end

function self = list_phase_species(self, LS)
    for ind=1:length(LS)
        Species = LS{ind};
        if ~self.strThProp.(Species).swtCondensed
           self.S.ind_nswt = [self.S.ind_nswt, ind];
        else
           self.S.ind_swt = [self.S.ind_swt, ind];
        end
    end
    % Number of gaseous species
    self.S.NG = length(self.S.ind_nswt);
end

function self = rearrange_species(self)
    self.S.LS = self.S.LS([self.S.ind_nswt, self.S.ind_swt]);
    if any(contains(self.S.LS,'Cbgrb'))
       self.S.ind_Cgr = find_ind({'Cbgrb'}, self.S.LS);
    end
    self.S.ind_nswt = [];
    self.S.ind_swt = [];
    % Index gaseous and condensed species
    self = list_phase_species(self, self.S.LS);
end

function self = Index_fixed_Species(self)
    % List of fixed gaseous species and Cgr
    self.S.ind_CO2 = find_ind({'CO2'}, self.S.LS);
    self.S.ind_CO  = find_ind({'CO'}, self.S.LS);
    self.S.ind_H2O = find_ind({'H2O'}, self.S.LS);
    self.S.ind_H2  = find_ind({'H2'}, self.S.LS);
    self.S.ind_O2  = find_ind({'O2'}, self.S.LS);
    self.S.ind_N2  = find_ind({'N2'}, self.S.LS);
    self.S.ind_He  = find_ind({'He'}, self.S.LS);
    self.S.ind_Ar  = find_ind({'Ar'}, self.S.LS);
    self.S.ind_Cgr  = find_ind({'Cbgrb'}, self.S.LS);
    self.S.ind_fixed = [self.S.ind_CO2,self.S.ind_CO,self.S.ind_H2O,...
        self.S.ind_H2,self.S.ind_O2,self.S.ind_N2,self.S.ind_He,...
        self.S.ind_Ar,self.S.ind_Cgr];
end

function self = Stoich_Matrix(self)
    self.C.A0.value = zeros(self.S.NS,self.E.NE);
    self.C.M0.value = zeros(self.S.NS,12);
    for i=1:self.S.NS
        txFormula = self.strThProp.(self.S.LS{i}).txFormula;
        self.strThProp.(self.S.LS{i}).Element_matrix = set_element_matrix(txFormula,self.E.elements);
        self.C.A0.value(i,self.strThProp.(self.S.LS{i}).Element_matrix(1,:)) = self.strThProp.(self.S.LS{i}).Element_matrix(2,:);
        self.C.M0.value(i,10) = self.strThProp.(self.S.LS{i}).swtCondensed;    
    end
    self.C.N0.value = self.C.M0.value(:, [1, 10]);
end

function self = Compute_minors_species(self)
    if self.Misc.FLAG_FIRST
        self.M.Lminors = length(self.M.minors_products);
        if self.M.Lminors > 0
            for n=self.M.Lminors:-1:1
                % Properties of other minor species under consideration, which can
                % be written in the generic form C_alpha H_beta O_gamma N_omega
                % Find index minor species
                self.M.ind_minor(n) = find(strcmp(self.S.LS, self.M.minors_products{n}));
            end
            self.C.alpha = self.C.A0.value(self.M.ind_minor, self.E.ind_C)';
            self.C.beta  = self.C.A0.value(self.M.ind_minor, self.E.ind_H)';
            self.C.gamma = self.C.A0.value(self.M.ind_minor, self.E.ind_O)';
            self.C.omega = self.C.A0.value(self.M.ind_minor, self.E.ind_N)';

            self.S.ind_all = [self.S.ind_fixed, self.M.ind_minor];
        else
            self.S.ind_all = self.S.ind_fixed;
        end
    end
end

function self = CH4_major(self)
    if any(contains(self.M.minors_products,'CH4')) && strcmp(self.PD.method, 'SEGREGATED')
        self.M.major_CH4 = true;
        self.M.ind_m_CH4 = find_ind({'CH4'}, self.M.minors_products);
    %     self.M.ind_m_C2H2= find_ind({'C2H2_acetylene'}, self.M.minors_products);
    %     self.M.ind_m_C6H6= find_ind({'C6H6'}, self.M.minors_products);

        self.M.ind_m_CH3 = find_ind({'CH3'}, self.M.minors_products);
        self.M.ind_m_H = find_ind({'H'}, self.M.minors_products);
        self.M.ind_m_CH = find_ind({'CH'}, self.M.minors_products);
    %     self.M.ind_m_C = find_ind({'C'}, self.M.minors_products);
    else
        self.M.major_CH4 = false;
    end
end

function self = OH_major(self) 
    if any(contains(self.M.minors_products,'OH')) && strcmp(self.PD.method, 'SEGREGATED')
        self.M.major_OH = true;
        self.M.ind_m_OH = find_ind({'OH'}, self.M.minors_products);
    else
        self.M.major_OH = false;
    end
end

function PT = Ask_problem(self)
    if self.Misc.FLAG_FOI
        fn = {'TP','HP','SP','TV','EV','SV','SHOCK_I','SHOCK_R','DET','DET_OVERDRIVEN'};
        [indx,tf] = listdlg('PromptString','Select a problem:',...
                                   'SelectionMode','single',...
                                   'ListString',fn,'ListSize',[150,120]);
        if tf
            PT = fn{indx};
        else
            error('Problem type not selected.')
        end
    else
        PT = self.PD.ProblemType;
    end
end