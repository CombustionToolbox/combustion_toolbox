T = linspace(1000, 3500, 5000); % [K]
P = 1;    % [atm]
% Combustion Toolbox
t_CT = test_CT(T, P);
% Cantera
t_CANTERA = test_CANTERA(T, P);
% test_SD(T,P)

% Print elapsed time
fprintf('\n COMPUTATIONAL COST\n');
fprintf('Time CT:      %1.4e seconds\n', t_CT);
fprintf('Time CANTERA: %1.4e seconds\n', t_CANTERA);
fprintf('Speed factor: %.2f times slower\n', t_CT / t_CANTERA);

%% SUB-PASS FUNCTIONS
function tEnd = test_CT(T,P)
    %% INITIALIZE
%     self = App({'H2','H','O','O2','OH','H2O','HO2','H2O2','C','CH','CH2','CH3','CH4','CO','CO2','HCO','CH2OH','CH3O','CH3OH','C2H','C2H4','C2H5','C2H6','HCCO','N','NH','NH2','NH3','NO','NO2','N2O','HNO','CN','HCN','HNCO','NCO','N2','C3H8'});
    self = App({'H2', 'H', 'O', 'O2', 'OH', 'H2O', 'HO2', 'H2O2'});
    %% INITIAL CONDITIONS
    self = set_prop(self, 'TR', 300, 'pR', 1 * 1.01325);
    self.PD.S_Fuel     = {'H2'};
    self.PD.N_Fuel     = 2;
    self.PD.S_Oxidizer = {'O2'};
%     self.PD.S_Inert    = {'N2'};
%     self.PD.proportion_inerts_O2 = 79/21;
    %% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
    self = set_prop(self, 'pP', P * 1.01325, 'TP', T); 
    %% TUNNING PROPERTIES
    self.TN.tolN = 1e-14;
    self.TN.tol0 = 1e-8;
    %% CONSTANTS (OPTIONAL)
    self.C.mintol_display = 1e-14;
    self.Misc.FLAG_RESULTS = false;
    %% SOLVE PROBLEM
    tStart = tic;
    self = solve_problem(self, 'TP');
    tEnd = toc(tStart);
    %% DISPLAY RESULTS (PLOTS)
    post_results(self);
end

function tEnd = test_CANTERA(T,P)
    mech = 'h2o2_bishnu1997.yaml';
    gas = Solution(mech);
    x = 'H2:2 O2:1';
    tStart = tic;
    for i = 1:length(T)
        set(gas,'Temperature',T(i),'Pressure',P * 101325.0,'MoleFractions',x);
        equilibrate(gas,'TP',0,1e-8,300,300);
    end
    tEnd = toc(tStart);
    xeq = moleFractions(gas);
    table(speciesNames(gas)', xeq);
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