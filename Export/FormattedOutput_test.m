function ExcelOutputMatrix = FormattedOutput_test(ExcelOutputMatrix,phi,str,NameSpecies)
R0 = 8.3144598;
% FormattedOutput_test(ExcelOutputMatrix,phi,T,p,P,hfR,DhTR,MassOrMolar,HeaddingOrData,strThProp)
% T must be specified in K
% p must be specified in bar
% MassOrMolar may be either 'mass' or 'molar'
% HeaddingOrData may be 'heading' or 'data'

% if nargin < 9
%     disp(['Error - Not enough arguments in formatted_output.m'])
% end

% set_elements;

% Ni = P(:,NE+1); % number of moles of species i in the mixture

% NM = sum(Ni);   % overall number of moles in the mixture
% Xi = Ni/NM;     % mass fraction of species i

% MM = 0;         % overall mass of mixture
% WM = 0;         % averge moleular mass of mixture
% NAS = 0;        % number of active species

% Yi_times_WM = zeros(size(Xi));

% fnm = fieldnames(strThProp);
% for i = 1:numel(fnm)
%     Wi = strThProp.(fnm{i}).mm/1000;
%     MM = MM + Wi*Ni(i);
%     WM = WM + Wi*Xi(i);
%     Yi_times_WM(i) = Wi*Xi(i);
%     if (Ni(i) > 0)|(any(strcmp(fnm{i},{'CO2','CO','H2O','H2','O2','N2'}))), NAS = NAS + 1; end
% end

% Yi = Yi_times_WM/WM;
% rho = p*101325/(Rg*T);

% hfP       = sum(P(:,NE+2));
% DhTP      = sum(P(:,NE+3));

% hM_molar  = (hfP + DhTP)/NM;
% hM_mass   = (hfP + DhTP)/MM;

% Dhf_molar = (hfP - hfR)/NM;
% Dhf_mass  = (hfP - hfR)/MM;

% Dh_molar = (hfP + DhTP - hfR - DhTR)/NM;
% Dh_mass  = (hfP + DhTP - hfR - DhTR)/MM;

% cPM_molar = sum(P(:,NE+6))/NM;
% cPM_mass  = sum(P(:,NE+6))/MM;

% cVM_molar = sum(P(:,NE+7))/NM;
% cVM_mass  = sum(P(:,NE+7))/MM;

% cVM_molar = cPM_molar - 8.314;
% cVM_mass  = cPM_mass - Rg;
%
% gamma     = cPM_molar/cVM_molar;
% gamma     = cPM_mass/cVM_mass;

% generation of screen/file output

% if strcmpi(MassOrMolar,'MASS')
%     composition_variable = 'Y';
composition_unit = 'kg';
% elseif strcmpi(MassOrMolar,'MOLAR')
composition_variable = 'X';
%     composition_unit = 'mol';
% else
%     disp(['Error - ',MassOrMolar,' unknown option for MassOrMolar input'])
% end

% if strcmpi(HeaddingOrData,'HEADING')

% heading output

ExcelOutputMatrix{1,1} = 'phi';
ExcelOutputMatrix{2,1} = '[-]';
ExcelOutputMatrix{1,2} = 'T';
ExcelOutputMatrix{2,2} = '[K]';
ExcelOutputMatrix{1,3} = 'p';
ExcelOutputMatrix{2,3} = '[bar]';
ExcelOutputMatrix{1,4} = 'W';
ExcelOutputMatrix{2,4} = '[g/mol]';
ExcelOutputMatrix{1,5} = 'Rg';
ExcelOutputMatrix{2,5} = ['[J/(',composition_unit,'·K)]'];
ExcelOutputMatrix{1,6} = 'rho';
ExcelOutputMatrix{2,6} = '[kg/m^3]';
ExcelOutputMatrix{1,7} = 'h';
ExcelOutputMatrix{2,7} = ['[kJ/',composition_unit,']'];
ExcelOutputMatrix{1,8} = 'DhT';
ExcelOutputMatrix{2,8} = ['[kJ/',composition_unit,']'];
ExcelOutputMatrix{1,9} = 'cP';
ExcelOutputMatrix{2,9} = ['[kJ/(',composition_unit,'·K)]'];
ExcelOutputMatrix{1,10} = 'gamma';
ExcelOutputMatrix{2,10} = '[-]';
ExcelOutputMatrix{1,11} = 's';
ExcelOutputMatrix{2,11} = ['[kJ/(',composition_unit,'·K)]'];
ExcelOutputMatrix{1,12} = 'u';
ExcelOutputMatrix{2,12} = '[m/s]';

SpeciesCounter = 0;
for j = 1:numel(NameSpecies)
    if (str{1}.Xi(j) > 0)||(any(strcmp(NameSpecies{j},{'CO2','CO','H2O','H2','O2','N2'})))
        SpeciesCounter = SpeciesCounter + 1;
        ExcelOutputMatrix{1,12+SpeciesCounter} = [composition_variable,NameSpecies{j}];
        ExcelOutputMatrix{2,12+SpeciesCounter} = '[-]';
    end
end

% ExcelOutputMatrix = [ExcelOutputMatrix; ExcelOutputMatrix];

% formatted_heading_1 = [''''];
% formatted_heading_2 = [''''];
% for i = 1:numel(ExcelOutputMatrix(1,:))
%     formatted_heading_1 = [formatted_heading_1, '%11s '];
%     formatted_heading_2 = [formatted_heading_2, '%11s '];
% end
% formatted_heading_1 = [formatted_heading_1, '\n'''];
% formatted_heading_2 = [formatted_heading_2, '\n'''];
% for i = 1:numel(ExcelOutputMatrix(1,:))
%     formatted_heading_1 = [formatted_heading_1, [',''',ExcelOutputMatrix{1,i},'''']];
%     formatted_heading_2 = [formatted_heading_2, [',''',ExcelOutputMatrix{2,i},'''']];
% end

%disp(formatted_heading_1)
%disp(formatted_heading_2)

% eval(['fprintf(',formatted_heading_1,')']);
% eval(['fprintf(',formatted_heading_2,')']);

% elseif strcmpi(HeaddingOrData,'DATA')

%% DATA OUTPUT
for i=length(phi):-1:1
    ExcelOutputMatrix{i+2,1} = phi(i);
    ExcelOutputMatrix{i+2,2} = str{i}.T;
    ExcelOutputMatrix{i+2,3} = str{i}.p;
    ExcelOutputMatrix{i+2,4} = str{i}.W;
    ExcelOutputMatrix{i+2,5} = R0/str{i}.W/1e-3;
    ExcelOutputMatrix{i+2,6} = str{i}.rho;
    ExcelOutputMatrix{i+2,7} = str{i}.h/str{i}.mi;
    ExcelOutputMatrix{i+2,8} = str{i}.DhT/str{i}.mi;
    ExcelOutputMatrix{i+2,9} = str{i}.cP/str{i}.mi*1e-3;
    %     eval(['ExcelOutputMatrixRow{1,7} = hM_',MassOrMolar,';']);
    %     eval(['ExcelOutputMatrixRow{1,8} = Dhf_',MassOrMolar,';']);
    %     eval(['ExcelOutputMatrixRow{1,9} = cPM_',MassOrMolar,';']);
    ExcelOutputMatrix{i+2,10} = str{i}.gamma;
    ExcelOutputMatrix{i+2,11} = str{i}.S;
    if isfield(str{i},'u')
        ExcelOutputMatrix{i+2,12} = str{i}.u;
    end
    SpeciesCounter = 0;
    for j = 1:numel(NameSpecies)
        if (str{i}.Xi(j) > 0)||(any(strcmp(NameSpecies{j},{'CO2','CO','H2O','H2','O2','N2'})))
            SpeciesCounter = SpeciesCounter + 1;
%             eval(['ExcelOutputMatrixRow{1,10+SpeciesCounter} = ',composition_variable,'i(',num2str(j),');']);
%             eval(['ExcelOutputMatrixRow{1,10+SpeciesCounter} = ',str{i}.Xi(j),';']);
            ExcelOutputMatrix{i+2,12+SpeciesCounter} = str{i}.Xi(j);
        end
    end
%     ExcelOutputMatrix = [ExcelOutputMatrix; ExcelOutputMatrix];
    
%     formatted_heading_1 = [''''];
%     for j = 1:numel(ExcelOutputMatrix(i,:))
%         formatted_heading_1 = [formatted_heading_1, '%11.3f '];
%     end
%     formatted_heading_1 = [formatted_heading_1, '\n'''];
%     for j = 1:numel(ExcelOutputMatrix(1,:))
%         formatted_heading_1 = [formatted_heading_1, [',',num2str(ExcelOutputMatrix{i+2,j})]];
%     end
    
    %disp(formatted_heading_1)
    
%     eval(['fprintf(',formatted_heading_1,')']);
end
%     else 
%     disp(['Error - ',HeaddingOrData,' unknown option for HeaddingOrData input']) 
% end