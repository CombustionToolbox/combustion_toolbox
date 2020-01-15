function [app,strThProp,strMaster] = Initialize()
%% LOAD/CALCULATE TABULATED DATA AND CONSTANTS
tic;
clearvars -except strMaster strThProp;
app.Misc.timer_0 = tic;
% path(path,pwd);
% path(path,[pwd,'\NASA_database']);
addpath(genpath(pwd))
% global factor_c Elements NE ind_C ind_H ind_O ind_N ind_He ind_Ar NameSpecies...
%     NSpecies displayRP displayT MassOrMolar filename mintol R0 strThProp...
%     minor_products CompleteOrIncomplete firstrow ProblemType mintol_display...
%     ExcelOutputMatrix save_Excel tolN idx_CO2 idx_CO idx_Cgr idx_H2O ...
%     idx_H2 idx_O2 idx_N2 idx_He idx_Ar FLAG_FIRST List_fixed_Species...
%     idx_fixed display_species A0
%% PROPERTIES DESCRIPTION
app.E.Description    = "Data of the chemical elements";
app.S.Description    = "Data of the chemical species";
app.M.Description    = "Data of minors products";
app.C.Description    = "Constants and tolerances";
app.C.A0.Description = "Stoichiometric matrix";
app.Misc.Description = "Miscelaneous properties";
app.PD.Description   = "Problem description";
app.PS.Description   = "Problem solution";
app.TN.Description   = "Tunning properties";
% Problem Description (PD):
app.PD.pR.Description  = "Temperature of the reactive species [K]";
app.PD.pR.Description  = "Pressure of the reactive species [bar]";
app.PD.phi.Description = "Equivalence ratio [-]";
%%
[app.E.Elements,~] = set_elements();
%% LOAD/CALCULATE TABULATED DATA
if ~exist('strMaster','var')
   if exist('strMaster.mat','file')
    disp('NASA database loaded from main path')
    load('strMaster_reduced.mat');
   else
       strMaster = ParseThermoInp_test(); % struct with tabulated data of all species.
       strMaster = strMaster_reduced(strMaster);
   end
else
    disp('NASA database already loaded')
end
%% LOAD/CALCULATE SHORT TABULATED DATA
if ~exist('strThProp','var')
   if exist('strThProp.mat','file')
    disp('NASA short database loaded from main path')
    load('strThProp.mat');
%     load('strThProp_HC_47.mat');
   else
       strThProp = GenerateShortDatabase_test(strMaster); % struct with tabulated data of selected species. 
   end
else
    disp('NASA short database already loaded')
end
%% MISCELANEOUS
app.Misc.FLAG_FIRST = true;
app.PD.Fuel.FLAG_phic = true;
%% INITIALIZE VALUES
app.M.display_species = {};
app.M.minor_products = {};
%% CONTAINED ELEMENTS
app.S.NameSpecies = fieldnames(strThProp); 
app.S.NSpecies = numel(app.S.NameSpecies);
% Match=cellfun(@(x) ismember(x,NameSpecies), Elements, 'UniformOutput', 0);
% Elements = Elements(cell2mat(Match));
% %%%%%%%%%%%%%%%%%%%%%%%%%
% Tmp = cell(size(NameSpecies));
% for k = 1:length(Tmp)
%     Specie = NameSpecies{k}; 
%     Specie(Specie>='0' & Specie<='9') = ' ';
%     idx = find([(Specie>='A' & Specie<='Z'), true]);
%     lgt = diff(idx);
%     Tmp{k} = strtrim(mat2cell(Specie, 1, lgt));
% end
% Elements = unique(cat(2,Tmp{:}));
% %%%%%%%%%%%%%%%%%%%%%%%%%
% Tmp = cell(size(NameSpecies));
for k = length(app.S.NameSpecies):-1:1
    Species = app.S.NameSpecies{k};
    % Change uppercase 'L' to  lowercase 'l'
    Species(strfind(Species,'AL')+1)='l';
    Species(strfind(Species,'CL')+1)='l';
    Species(strfind(Species,'TL')+1)='l';
    Species(strfind(Species,'FL')+1)='l';
    % -----------------------------------------------
    Species(Species>='0' & Species<='9') = ' ';
    
    [idx0,idxf] = regexp(Species,"minus"); Species(idx0:idxf) = ' ';
    [idx0,idxf] = regexp(Species,"plus"); Species(idx0:idxf) = ' ';
    
    idx = find([(Species>='A' & Species<='Z'), true]);
    lgt = diff(idx);
    Tmp{k,1} = strtrim(mat2cell(Species, 1, lgt));
end
aux = unique(cat(2,Tmp{:}));
n_pass = [];
for n=length(aux):-1:1
    if any(strcmp(aux(n),app.E.Elements)) % Check Element existence
        n_pass = [n_pass, n];
    end
end   
app.E.Elements = aux(n_pass);
app.E.NE = numel(app.E.Elements);
%% INDEX OF EVALUABLE ELEMENTS
app.E.ind_C = find(strcmp(app.E.Elements,'C'));
app.E.ind_H = find(strcmp(app.E.Elements,'H'));
app.E.ind_O = find(strcmp(app.E.Elements,'O'));
app.E.ind_N = find(strcmp(app.E.Elements,'N'));
app.E.ind_He = find(strcmp(app.E.Elements,'He'));
app.E.ind_Ar = find(strcmp(app.E.Elements,'Ar'));
%% INDEX OF EVALUABLE SPECIES
app.S.idx_CO2 = find_idx({'CO2'},app.S.NameSpecies);
app.S.idx_CO  = find_idx({'CO'},app.S.NameSpecies);
app.S.idx_H2O = find_idx({'H2O'},app.S.NameSpecies);
app.S.idx_H2  = find_idx({'H2'},app.S.NameSpecies);
app.S.idx_O2  = find_idx({'O2'},app.S.NameSpecies);
app.S.idx_N2  = find_idx({'N2'},app.S.NameSpecies);
app.S.idx_He  = find_idx({'He'},app.S.NameSpecies);
app.S.idx_Ar  = find_idx({'Ar'},app.S.NameSpecies);
app.S.idx_Cgr = find_idx({'Cbgrb'},app.S.NameSpecies);
app.S.List_fixed_Species = {'CO2','CO','H2O','H2','O2','N2','He','Ar','Cbgrb'};
app.S.idx_fixed = [app.S.idx_CO2,app.S.idx_CO,app.S.idx_H2O,app.S.idx_H2,app.S.idx_O2,app.S.idx_N2,app.S.idx_He,app.S.idx_Ar,app.S.idx_Cgr];
%% CONSTANTS & TOLERANCES
app.C.R0 = 8.3144598; % [J/(K mol)]. Universal gas constant
app.C.MassOrMolar = 'mass';
app.C.filename = 'output';
app.C.firstrow = 1;
app.C.ExcelOutputMatrix = [];
app.C.mintol_display = 1e-8;
app.C.mintol = 1e-5;
app.C.tolN = 1e-14;
app.C.FLAG_FIRST = true;
%% Element_matrix
app.C.A0.Value = zeros(app.S.NSpecies,app.E.NE);
app.C.M0.Value = zeros(app.S.NSpecies,12);
for i=1:app.S.NSpecies
    txFormula = strThProp.(app.S.NameSpecies{i,1}).txFormula;
    strThProp.(app.S.NameSpecies{i,1}).Element_matrix = set_element_matrix(txFormula,app.E.Elements);
    app.C.A0.Value(i,strThProp.(app.S.NameSpecies{i,1}).Element_matrix(1,:)) = strThProp.(app.S.NameSpecies{i,1}).Element_matrix(2,:);
    app.C.M0.Value(i,10) = strThProp.(app.S.NameSpecies{i,1}).swtCondensed;    
end
%% GUESS INITIAL CALCULATION
app.TN.guess = [2000,2000,0,1.5,2];
% TOLERANCES
app.TN.ERRFT = 1e-5;
app.TN.ERRFU = 1e-5;
toc;