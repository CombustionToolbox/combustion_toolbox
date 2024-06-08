% -------------------------------------------------------------------------
% EXAMPLE: ROCKET Propellants considering a Finite-Area-Chamber (FAC)
%
% Compute adiabatic temperature and equilibrium composition at constant
% pressure (e.g., 1.01325 bar) for lean to rich LH2-LOX mixtures at
% standard conditions, a set of 24 species considered and a set of
% equivalence ratios phi contained in (1, 5) [-]
%   
% HYDROGEN_L == {'H','H2O','OH','H2','O','O3','O2','HO2','H2O2',...
%                'H2bLb','O2bLb'}
%   
% See wiki or setListspecies method from ChemicalSystem class for more
% predefined sets of species
%
% @author: Alberto Cuadra Lara
%          Postdoctoral researcher - Group Fluid Mechanics
%          Universidad Carlos III de Madrid
%                 
% Last update July 22 2022
% -------------------------------------------------------------------------

FLAG_SAVE = true;
%% SOLVE PROBLEM
phi = 1:0.01:4;
Aratio = 1.1:0.1:2.6;
DisplaySpecies = {'CO2','CO','H2O','H2','O2','C2H2_acetylene',...
          'C2H4','C2H6','CH2CO_ketene','CH3','CH3CHO_ethanal','CH3OH',...
          'CH4','COOH','H','H2O2','HCHO_formaldehy','HCO','HCOOH','HO2',...
          'O','OH','Cbgrb'};

for i = length(Aratio):-1:1
    %% INITIALIZE
    self = App('HC/O2/N2 PROPELLANTS');
    self.TN.tolN = 1e-18;
    self.Misc.display_species = DisplaySpecies;
    %% INITIAL CONDITIONS
    self.PD.S_Fuel     = {'RP_1'};
    self.PD.S_Oxidizer = {'O2bLb'};
    self.PD.FLAG_IAC = false;
    %% ADDITIONAL INPUTS (DEPENDS OF THE PROBLEM SELECTED)
    self = set_prop(self, 'Aratio_c', 2, 'Aratio', Aratio(i));
    self = set_prop(self, 'TR', 298.15, 'pR', 22, 'phi', phi);
    self = solve_problem(self, 'ROCKET');
    Z(:, i) = cell2vector(self.PS.strP, 'I_sp');
    OF(:, i) = cell2vector(self.PS.strR, 'OF');
    %% DISPLAY RESULTS (PLOTS)
    if FLAG_SAVE
        % Plot molar fraction
        self.C.mintol_display = 1e-16;
        self.Misc.config.labelx = 'Mixture ratio $O/F$';
        self.Misc.config.labely = 'Molar fraction $X_i$';
        self.Misc.config.title = strcat('Area ratio $A_{\rm exit}/A_{\rm throat} = ', sprintf('%.2f', Aratio(i)), '$');
        ax = displaysweepresults(self, self.PS.strP, OF(:, i));
        set_title(ax, self.Misc.config);
        % Save plots
        folderpath = strcat(pwd,'\Validations\Figures\');
        stack_trace = dbstack;
        filename = strcat('RP1_LOX_', sprintf('%d', i));
        saveas(gcf, strcat(folderpath, filename, '_molar'), 'svg');
    end
end
% Contour plot
[X, Y] = meshgrid(OF(:, i), Aratio);
ax = set_figure;
contourf(ax, X, Y, Z', 10);
c = colorbar;
xlabel(ax, 'Mixture ratio $O/F$');
ylabel(ax, 'Area ratio $A_{\rm exit}/A_{\rm throat}$');
c.Label.String = 'Specific impulse at sea level $I_{sp}$';
c.Label.Interpreter = 'latex';
set(ax, 'Layer','top')