function gui_update_UITable_R(obj, app)
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

% SUB-PASS FUNCTIONS
function typeSpecies = get_typeSpecies(app)
    % Create cell array with the type of species in the mixture
    typeFuel     = create_cell_ntimes('Fuel', length(app.PD.N_Fuel));
    typeOxidizer = create_cell_ntimes('Oxidizer', length(app.PD.N_Oxidizer));
    typeInert    = create_cell_ntimes('Inert', length(app.PD.N_Inert));
    typeSpecies  = [typeInert, typeOxidizer, typeFuel];
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