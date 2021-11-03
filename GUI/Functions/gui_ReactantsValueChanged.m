function gui_ReactantsValueChanged(obj, event)
    % Update values of the UITable items:
    % 1. with a given predefined set of reactants
    % 2. with the new species added in the finder
    try
        % Initialize app (fast: transfer DB)
        app = App('fast', obj.DB_master, obj.DB);
        % Update reactants & GUI
        obj = gui_update_Reactants(obj, event, app);
        % Update equivalence ratio 
        obj.edit_phi2.Value = obj.edit_phi.Value;
    catch ME
      errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
      ME.stack(1).name, ME.stack(1).line, ME.message);
      fprintf('%s\n', errorMessage);
      uiwait(warndlg(errorMessage));
    end
end

% SUB-PASS FUNCTIONS
function obj = gui_update_Reactants(obj, event, app)
    FLAG_IDEAL_AIR = obj.IdealAirCheckBox.Value;
    if strcmp(obj.Reactants.Value, '1') % No species selected
        gui_empty_Reactants(obj);
        return
    end
    % Default value of equivalence ratio is 1
    % Set species in the mixture
    switch obj.Reactants.Value
        case '2' % Air
            app = set_air(app, FLAG_IDEAL_AIR);
        case '3' % Methane + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'CH4'};
        case '4' % Ethane + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'C2H6'};
        case '5' % Propane + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'C3H8'};
        case '6' % Acetylene + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'C2H2_acetylene'};
        case '7' % Ethylene + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'C2H4'};
        case '8' % Bencene + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'C6H6'};
        case '9' % Iso-octane + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'C8H18_isooctane'};
        case '10' % Hydrogen + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'H2'};
        case '11' % Carbon monoxide + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'CO'};
        case '12' % Methanol + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'CH3OH'};
        case '13' % Ethanol + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'C2H5OH'};
        case '14' % Natural Gas + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'CH4','C2H6','C3H8'};
            app.PD.N_Fuel = [0.85, 0.1, 0.05];
        case '15' % Syngas + Air
            app = set_air(app, FLAG_IDEAL_AIR);
            app.PD.S_Fuel = {'CO','H2'};  
            app.PD.N_Fuel = [0.5, 0.5];
        otherwise % SET NEW SPECIES
            try
                species = gui_seeker_species(obj, event);
            catch
                message = {'Species not found.'};
                uialert(obj.UIFigure, message, 'Warning', 'Icon', 'warning');
                return
            end
            % Get data of the current mixture
            if ~isempty(obj.UITable_R.Data)
                app = get_current_reactants_gui(obj, app, current_species, current_moles);
            end
            % Add new species to the mixture (fuel by default)
            app.PD.S_Fuel = [app.PD.S_Fuel, {species}];
            app.PD.N_Fuel = [app.PD.N_Fuel, 0];
    end
    % Compute properties of the mixture
    [obj, app] = gui_compute_propReactants(obj, app);
    % Update UITable classes
    obj = gui_update_UITable_R(obj, app);
    obj.UITable_P.Data = obj.UITable_R.Data(:, 1);    % (species, numer of moles, mole fractions, temperature)
    obj.UITable_R2.Data = obj.UITable_R.Data(:, 1:3); % (species, numer of moles, mole fractions)
end

% SUB-PASS FUNCTIONS
function obj = gui_update_UITable_R(obj, app)
    % Update data in the UITable_R with the next order: Inert -> Oxidizer -> Fuel
    species = [app.PD.S_Inert, app.PD.S_Oxidizer, app.PD.S_Fuel];
    Nspecies = length(species); 
    if isempty(app.PD.S_Fuel)
        app.PD.N_Fuel = []; % Set to 1 by default
    end
    moles = [app.PD.N_Inert, app.PD.N_Oxidizer, app.PD.N_Fuel];
    molar_fraction = moles/sum(moles); % It is easier to recompute
    typeSpecies = get_typeSpecies(app);
    temperature =  create_cell_ntimes(app.PD.TR.value, Nspecies);
    obj.UITable_R.Data = [species; vector2cell(moles); vector2cell(molar_fraction); typeSpecies; temperature]';
end

function C = vector2cell(value)
    % Create cell array from vector
    for i=length(value):-1:1
        C(i) = {value(i)};
    end
end

function C = create_cell_ntimes(varargin)
    % Create cell array with the same string n-times
    if nargin > 2
        value = varargin{1};
        C = varargin{3};
    elseif nargin > 1
        value = varargin{1};
        n = varargin{2};
        C = cell(1, n);
    else
        error('Error sub-pass fuinction @create_cell_ntimes inside @guiReactantsValueChanged');
    end
    % Set value
    C(:) = {value};
end

function typeSpecies = get_typeSpecies(app)
    % Create cell array with the type of species in the mixture
    typeFuel     = create_cell_ntimes('Fuel', length(app.PD.N_Fuel));
    typeOxidizer = create_cell_ntimes('Oxidizer', length(app.PD.N_Oxidizer));
    typeInert    = create_cell_ntimes('Inert', length(app.PD.N_Inert));
    typeSpecies  = [typeInert, typeOxidizer, typeFuel];
end

function app = set_air(app, FLAG_IDEAL_AIR)
    % Incluide air in the mixture
    if FLAG_IDEAL_AIR
        app.PD.S_Oxidizer = {'O2'};
        app.PD.S_Inert = {'N2'};
        app.PD.proportion_inerts_O2 = 79/21;
    else
        app.PD.S_Oxidizer = {'O2'};
        app.PD.S_Inert = {'N2', 'Ar', 'CO2'};
        app.PD.proportion_inerts_O2 = [78.084, 0.9365, 0.0319] ./ 20.9476;
    end
end

function obj = gui_empty_Reactants(obj)
    % Clear data UITables and set to default the value of the equivalence ratio (-)
    obj.UITable_R.Data  = [];
    obj.UITable_P.Data  = [];
    obj.UITable_R2.Data = [];
    obj.edit_phi.Value = '-';
end