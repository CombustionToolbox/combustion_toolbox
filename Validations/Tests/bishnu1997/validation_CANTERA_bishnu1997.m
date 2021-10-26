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
    self = App({'H2','H2O','H2O2','H','OH','HO2','O2','O'});
    %% INITIAL CONDITIONS
    self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325, 'phi', 1);
    self.PD.S_Fuel     = {'H2'};
    self.PD.S_Oxidizer = {'O2'};
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
    mech = 'h2o2_bishnu1997.yaml';
    gas = Solution(mech);
    x = 'H2:2 O2:1';
    set(gas,'Temperature',T,'Pressure',P * 101325.0,'MoleFractions',x);
    tStart = tic;
    cj_speed = CJspeed(P, T, x, mech);
    gas = PostShock_eq(cj_speed,P, T, x, mech)
    tEnd = toc(tStart);
    xeq = moleFractions(gas);
    table(speciesNames(gas)', xeq)
end