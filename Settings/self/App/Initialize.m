function self = Initialize(self)
    % Check if minors products species are contained in DB
    [self.strThProp, self.E, self.S, self.C] = check_database(self, self.strMaster, self.strThProp);
    % Sort species: first gaseous species, secondly condensed species
    self = list_phase_species(self, self.S.LS);
    % Stoichiometric Matrix
    self = Stoich_Matrix(self);
    % Ask problem type
    self.PD.ProblemType = Ask_problem(self);
end


%% NESTED FUNCTIONS
function self = list_phase_species(self, LS)
    for ind=1:length(LS)
        Species = LS{ind};
        if ~self.strThProp.(Species).swtCondensed
           self.S.ind_nswt = [self.S.ind_nswt, ind];
        else
           self.S.ind_swt = [self.S.ind_swt, ind];
           if ~self.strThProp.(Species).ctTInt
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
        txFormula = self.strThProp.(self.S.LS{i}).txFormula;
        self.strThProp.(self.S.LS{i}).Element_matrix = set_element_matrix(txFormula,self.E.elements);
        self.C.A0.value(i,self.strThProp.(self.S.LS{i}).Element_matrix(1,:)) = self.strThProp.(self.S.LS{i}).Element_matrix(2,:);
        self.C.M0.value(i,10) = self.strThProp.(self.S.LS{i}).swtCondensed;    
    end
    self.C.N0.value = self.C.M0.value(:, [1, 10]);
end


function PT = Ask_problem(self)
    if self.Misc.FLAG_FOI && ~self.Misc.FLAG_GUI
        try
        fn = {'TP','HP','SP','TV','EV','SV','SHOCK_I','SHOCK_R','DET','DET_OVERDRIVEN'};
        [indx, ~] = listdlg('PromptString','Select a problem:',...
                                   'SelectionMode','single',...
                                   'ListString',fn,'ListSize',[150,150]);
        PT = fn{indx};
        catch
            error('Problem type not selected.')
        end
    else
        if self.Misc.FLAG_GUI
            self = get_PT_GUI(self);
        end
        PT = self.PD.ProblemType;
    end
end

function self = get_PT_GUI(self)
    switch self.ProblemType.Value
        case '1' % * TP: Equilibrium composition at defined T and p
            self.PD.ProblemType = 'TP';
        case '2' % * HP: Adiabatic T and composition at constant p
            self.PD.ProblemType = 'HP';
        case '3' % * SP: Isentropic (i.e., adiabatic) compression/expansion to a specified p
            self.PD.ProblemType = 'SP';
        case '4' % * TV: Equilibrium composition at defined T and constant v
            self.PD.ProblemType = 'TV';  
        case '5' % * EV: Equilibrium composition at Adiabatic T and constant v
            self.PD.ProblemType = 'EV';
        case '6' % * SV: Isentropic (i.e., fast adiabatic) compression/expansion to a specified v
            self.PD.ProblemType = 'SV';
        case '7' % * SHOCK_I: CALCULATE PLANAR INCIDENT SHOCK WAVE
            self.PD.ProblemType = 'SHOCK_I';
        case '8' % * SHOCK_R: CALCULATE PLANAR POST-REFLECTED SHOCK STATE
            self.PD.ProblemType = 'SHOCK_R';
        case '9' % * DET: CALCULATE CHAPMAN-JOUGET STATE (CJ UPPER STATE)
            self.PD.ProblemType = 'DET';
        case '10' % * DET: CALCULATE OVERDRIVEN CHAPMAN-JOUGET STATE (CJ UPPER STATE)
            self.PD.ProblemType = 'DET_OVERDRIVEN';
    end
end