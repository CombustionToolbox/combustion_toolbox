T = 300; % [K]
P = 1;    % [atm]
% Combustion Toolbox
t_CT = test_CT(T, P);
% Cantera
% t_CANTERA = test_CANTERA(T, P);
test_SD(T,P)

%% SUB-PASS FUNCTIONS
function tEnd = test_CT(T,P)
    %% INITIALIZE
    self = App({'H2','H','O','O2','OH','H2O','HO2','H2O2','C','CH','CH2','CH3','CH4','CO','CO2','HCO','CH2OH','CH3O','CH3OH','C2H','C2H4','C2H5','C2H6','HCCO','N','NH','NH2','NH3','NO','NO2','N2O','HNO','CN','HCN','HNCO','NCO','N2','C3H8'});
    %% INITIAL CONDITIONS
    self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325, 'phi', 1);
    self.PD.S_Fuel     = {'CH4'};
    self.PD.S_Oxidizer = {'O2'};
    self.PD.S_Inert    = {'N2'};
    self.PD.proportion_inerts_O2 = 79/21;
    %% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
    self = set_prop(self, 'pP', P * 1.01325, 'TP', T); 
    %% TUNNING PROPERTIES
    self.TN.tolN = 1e-14;
    self.tol0 = 1e-8;
    %% CONSTANTS (OPTIONAL)
    self.C.mintol_display = 1e-14;
    %% SOLVE PROBLEM
    tStart = tic;
    self = SolveProblem(self, 'DET');
    tEnd = toc(tStart);
    %% DISPLAY RESULTS (PLOTS)
    postResults(self);
    
end

function tEnd = test_CANTERA(T,P)
    mech = 'h2o2_bishnu1997.yaml';
    gas = Solution(mech);
    x = 'H2:2 O2:1';
    set(gas,'Temperature',T,'Pressure',P * 101325.0,'MoleFractions',x);
    tStart = tic;
    equilibrate(gas,'HP',0,1e-8,300,300)
    tEnd = toc(tStart);
    xeq = moleFractions(gas);
    table(speciesNames(gas)', xeq)
end

function tEnd = test_SD(T,P)
    P = P * 1e5;
    mech = 'gri30.yaml';
    x = 'CH4:1 O2:2 N2:7.52';
    tStart = tic;
    cj_speed = CJspeed(P, T, x, mech);
    gas = PostShock_eq(cj_speed,P, T, x, mech)
    tEnd = toc(tStart);
    xeq = moleFractions(gas);
    table(speciesNames(gas)', xeq)
end