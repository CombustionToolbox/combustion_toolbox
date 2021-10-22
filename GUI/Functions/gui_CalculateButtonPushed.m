function self = gui_CalculateButtonPushed(app)
%     % Start
%     self.Lamp.Color = 'Yellow';
% %     self = self(self, self.UITable_P.Data(:, 1)');
%     self = App(self, 'Soot Formation');
%     % Get initial parameters
%     self = gui_get_parameters(self);
%     % Constants
%     self.C.l_phi = length(self.PD.phi.value);
%     self.ind_Fuel = strcmp(self.UITable_R.Data(:,4), 'Fuel');
%     self.ind_Oxidizer = strcmp(self.UITable_R.Data(:,4), 'Oxidant');
%     self.ind_Inert = strcmp(self.UITable_R.Data(:,4), 'Inert');
%     % Get Name Reactants species
%     gui_get_reactants(self);
%     for i=self.C.l_phi:-1:1 % Evading preallocate struct
%         if self.flag_PR1; self.PD.TR.value = self.PR1_vector(i); end
%         if self.flag_PR2; self.PD.pR.value = self.PR2_vector(i); end
%         if self.flag_PP1; self.PD.TP.value = self.PP1_vector(i); end
%         if self.flag_PP2; self.PD.pP.value = self.PP2_vector(i); end
%         % Define Fuel
%         self = Define_F(self);
%         % Define Oxidizer
%         self = Define_O(self);
%         % Define Inert
%         self = Define_I(self);
%         % Compute properties
%         self = Define_FOI(self, i);
%         % Solve Problem selected
%         self = SolveProblem(self, i);
%         % Display results command window
%         results(self, i);
%     end
%     self.Lamp.Color = 'Green';
% end

% SUB-PASS FUNCTIONS


self = App(app, 'Soot formation');
%% PROBLEM CONDITIONS
[self.PD.TR.value, flag_PR1] = gui_get_prop(app, app.PR1.Value, 'TR');
[self.PD.pR.value, flag_PR2] = gui_get_prop(app, app.PR2.Value, 'pR');
[self.PD.phi.value, flag_phi] = gui_get_prop(app, app.edit_phi.Value, 'phi');
%% PROBLEM TYPE
switch self.PD.ProblemType
    case 'TP' % * TP: Equilibrium composition at defined T and p
        [self.PD.TP.value, flag_PP1] = gui_get_prop(app, app.PP1.Value, 'TP');
        [self.PD.pP.value, flag_PP2] = gui_get_prop(app, app.PP2.Value, 'pP');
    case 'HP' % * HP: Adiabatic T and composition at constant p
        self.PD.pP.value = self.PD.pR.value;
    case 'SP' % * SP: Isentropic (i.e., adiabatic) compression/expansion to a specified p
        self.PD.pP.value = 10:1:50; self.PD.phi.value = 1*ones(1, length(self.PD.pP.value));
%         app.PD.pP.value = 10*ones(1, length(app.PD.phi.value));
    case 'TV' % * TV: Equilibrium composition at defined T and constant v
        self.PD.TP.value = 2000;
        self.PD.pP.value = self.PD.pR.value; % guess
    case 'EV' % * EV: Equilibrium composition at Adiabatic T and constant v
        self.PD.pP.value = self.PD.pR.value;
        % app.PD.pR.value = logspace(0,2,20); app.PD.phi.value = 1*ones(1,length(app.PD.pR.value));
    case 'SV' % * SV: Isentropic (i.e., fast adiabatic) compression/expansion to a specified v
        % REMARK! vP_vR > 1 --> expansion, vP_vR < 1 --> compression
        self.PD.vP_vR.value = 0.5:0.01:2; self.PD.phi.value = 1*ones(1, length(self.PD.vP_vR.value));
    case 'SHOCK_I' % * SHOCK_I: CALCULATE PLANAR INCIDENT SHOCK WAVE
        u1 = logspace(2, 5, 500);
        u1 = u1(u1<20000); u1 = u1(u1>=360);
%         u1 = [356,433,534,658,811,1000,1233,1520,1874,2310,2848,3511,4329,5337,6579,8111,9500,12328,15999,18421,21210,24421,28118,32375,37276,42919,49417,56899,65513];
%         u1 = linspace(360, 9000, 1000);
%         u1 = 20000;
        self.PD.u1.value = u1; self.PD.phi.value = ones(1,length(self.PD.u1.value));
    case 'SHOCK_R' % * SHOCK_R: CALCULATE PLANAR POST-REFLECTED SHOCK STATE
        u1 = linspace(400, 6000, 1000);
%         u1 = 2000;
        self.PD.u1.value = u1; self.PD.phi.value = ones(1,length(self.PD.u1.value));
    case 'DET' % * DET: CALCULATE CHAPMAN-JOUGET STATE (CJ UPPER STATE)
%         app.PD.TR_vector.value = app.PD.TR.value;
    case 'DET_OVERDRIVEN' % * DET_OVERDRIVEN: CALCULATE OVERDRIVEN DETONATION
        self.PD.overdriven.value  = 1:0.1:10; self.PD.phi.value = 1*ones(1,length(self.PD.overdriven.value));
end
%% LOOP
self.C.l_phi = length(self.PD.phi.value);
tic
for i=self.C.l_phi:-1:1
% DEFINE FUEL
self.PD.S_Fuel = {'CH4'}; self.PD.N_Fuel = 1;
self = Define_F(self);
% DEFINE OXIDIZER
self.PD.S_Oxidizer = {'O2'}; self.PD.N_Oxidizer = self.PD.phi_t/self.PD.phi.value(i);
self = Define_O(self);
% DEFINE DILUENTS/INERTS
self.PD.proportion_N2_O2 = 79/21;
self.PD.S_Inert = {'N2'}; self.PD.N_Inert = self.PD.phi_t/self.PD.phi.value(i) * self.PD.proportion_N2_O2;
self = Define_I(self);
% COMPUTE PROPERTIES
self = Define_FOI(self, i);
% PROBLEM TYPE
self = SolveProblem(self, i);
% DISPLAY RESULTS COMMAND WINDOW
results(self, i);
end
toc
%% DISPLAY RESULTS (PLOTS)
self.Misc.display_species = {};
% app.Misc.display_species = {'CO','CO2','H','HO2','H2','H2O','NO','NO2','N2','O','OH','O2','Cbgrb'};
closing(self);
end