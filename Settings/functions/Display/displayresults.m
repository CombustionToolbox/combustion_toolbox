function displayresults(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   varargin = (strR,str2,strP) or varargin = (strR,strP)
%   mix1  = Prop. state 1 (phi,species,...)
%   mix2  = Prop. state 2 (phi,species,...)
%   mix3  = Prop. state 3 (phi,species,...)
%   ...
% OUTPUT:
%   results on command window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% help displayresults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ProblemType = varargin{end-2};
mintol_display = varargin{end-1};
ListSpecies = varargin{end};

if nargin == 5
    mix1 = varargin{1}; mix2 = varargin{2};
    
    fprintf('***********************************************************\n');
    fprintf('-----------------------------------------------------------\n');
    fprintf('Problem type: %s  | phi = %4.4f\n',ProblemType, equivalenceRatio(mix1));
    fprintf('-----------------------------------------------------------\n');
    if contains(ProblemType, 'SHOCK') || contains(ProblemType, 'DET')
        fprintf('               |    STATE 1      |       STATE 2\n');
    else
        fprintf('               |    REACTANTS    |      PRODUCTS\n');
    end
    
    fprintf('T [K]          |   %12.4f  |   %12.4f\n', temperature(mix1), temperature(mix2));
    fprintf('p [bar]        |   %12.4f  |   %12.4f\n', pressure(mix1), pressure(mix2));
    fprintf('r [kg/m3]      |   %12.4f  |   %12.4f\n', density(mix1), density(mix2));
    fprintf('h [kJ/kg]      |   %12.4f  |   %12.4f\n', enthalpy_mass(mix1), enthalpy_mass(mix2));
    fprintf('e [kJ/kg]      |   %12.4f  |   %12.4f\n', intEnergy_mass(mix1), intEnergy_mass(mix2));
    fprintf('g [kJ/kg]      |   %12.4f  |   %12.4f\n', gibbs_mass(mix1), gibbs_mass(mix2));
    fprintf('s [kJ/(kg-K)]  |   %12.4f  |   %12.4f\n', entropy_mass(mix1), entropy_mass(mix2));
    fprintf('W [g/mol]      |   %12.4f  |   %12.4f\n', meanMolecularWeight(mix1), meanMolecularWeight(mix2));
    fprintf('(dlV/dlp)T [-] |                 |   %12.4f\n', mix2.dVdp_T);
    fprintf('(dlV/dlT)p [-] |                 |   %12.4f\n', mix2.dVdT_p);
    fprintf('cp [kJ/(kg-K)] |   %12.4f  |   %12.4f\n', cp_mass(mix1), cp_mass(mix2));
    fprintf('gamma [-]      |   %12.4f  |   %12.4f\n', adiabaticIndex(mix1), adiabaticIndex_sound(mix2));
    fprintf('sound vel [m/s]|   %12.4f  |   %12.4f\n', soundspeed(mix1), soundspeed(mix2));
    if contains(ProblemType, 'SHOCK') || contains(ProblemType,'DET')
        fprintf('u [m/s]        |   %12.4f  |   %12.4f\n', velocity_relative(mix1), mix2.v_shock);
        fprintf('Mach number [-]|   %12.4f  |   %12.4f\n', velocity_relative(mix1)/soundspeed(mix1), mix2.v_shock/soundspeed(mix2));
    end
    if contains(ProblemType, '_OBLIQUE') || contains(ProblemType, '_POLAR')
        fprintf('------------------------------------------------------------------------\n');
        fprintf('PARAMETERS\n');
        fprintf('min wave  [deg]|                 |   %12.4f\n', mix2.beta_min);
        if contains(ProblemType, '_OBLIQUE')
            fprintf('wave angle[deg]|                 |   %12.4f\n', mix2.beta);
            fprintf('deflection[deg]|                 |   %12.4f\n', mix2.theta);
        else
            fprintf('max def.  [deg]|                 |   %12.4f\n', mix2.theta_max);
            fprintf('sonic def.[deg]|                 |   %12.4f\n', mix2.theta_sonic);
        end
    end
    fprintf('-----------------------------------------------------------\n');
    fprintf('REACTANTS               Xi [-]\n');
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mix1.Xi(:), ind_sort] = sort(mix1.Xi(:), 'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = mix1.Xi>mintol_display;
    Xminor = sum(mix1.Xi(~j));
    for i=1:length(j)
        if j(i)
            fprintf('%-20s %1.4e\n', ListSpecies{ind_sort(i)}, mix1.Xi(i));
        end
    end
    Nminor = length(mix1.Xi)-sum(j);
    s_space_Nminor = char(32 * ones(1, 4 - numel(num2str(Nminor))));
    fprintf('MINORS[+%d] %s     %12.4e\n\n', Nminor, s_space_Nminor, Xminor);
    fprintf('TOTAL            %14.4e\n',sum(mix1.Xi));
    fprintf('-----------------------------------------------------------\n');
    fprintf('PRODUCTS                Xi [-]\n');
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mix2.Xi(:), ind_sort] = sort(mix2.Xi(:),'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = mix2.Xi>mintol_display;
    Xminor = sum(mix2.Xi(~j));
    for i=1:length(j)
        if j(i)
            fprintf('%-20s %1.4e\n', ListSpecies{ind_sort(i)}, mix2.Xi(i));
        end
    end
    Nminor = length(mix2.Xi) - sum(j);
    s_space_Nminor = char(32 * ones(1, 4 - numel(num2str(Nminor))));
    fprintf('MINORS[+%d] %s     %12.4e\n\n', Nminor, s_space_Nminor, Xminor);
    fprintf('TOTAL            %14.4e\n',sum(mix2.Xi));
    fprintf('------------------------------------------------------------------------\n');
    fprintf('************************************************************************\n\n\n');
elseif nargin == 6
    mix1 = varargin{1}; mix2 = varargin{2}; mix3 = varargin{3};
    fprintf('************************************************************************\n');
    fprintf('------------------------------------------------------------------------\n');
    fprintf('Problem type: %s  | phi = %4.4f\n',ProblemType, equivalenceRatio(mix1));
    fprintf('------------------------------------------------------------------------\n');
    if contains(ProblemType, '_R')
        fprintf('               |     STATE 1     |     STATE 2     |      STATE 3\n');
    elseif contains(ProblemType, '_OBLIQUE')
        fprintf('               |     STATE 1     |     STATE 2-W   |      STATE 2-S\n');
    else
        fprintf('               |  INLET CHAMBER  | OUTLET CHAMBER  |      THROAT \n');
    end
    
    fprintf('T [K]          |   %12.4f  |   %12.4f  |   %12.4f\n', temperature(mix1), temperature(mix2), temperature(mix3));
    fprintf('p [bar]        |   %12.4f  |   %12.4f  |   %12.4f\n', pressure(mix1), pressure(mix2), pressure(mix3));
    fprintf('r [kg/m3]      |   %12.4f  |   %12.4f  |   %12.4f\n', density(mix1), density(mix2), density(mix3));
    fprintf('h [kJ/kg]      |   %12.4f  |   %12.4f  |   %12.4f\n', enthalpy_mass(mix1), enthalpy_mass(mix2), enthalpy_mass(mix3));
    fprintf('e [kJ/kg]      |   %12.4f  |   %12.4f  |   %12.4f\n', intEnergy_mass(mix1), intEnergy_mass(mix2), intEnergy_mass(mix3));
    fprintf('g [kJ/kg]      |   %12.4f  |   %12.4f  |   %12.4f\n', gibbs_mass(mix1), gibbs_mass(mix2),  gibbs_mass(mix3));
    fprintf('s [kJ/(kg-K)]  |   %12.4f  |   %12.4f  |   %12.4f\n', entropy_mass(mix1), entropy_mass(mix2), entropy_mass(mix3));
    fprintf('W [g/mol]      |   %12.4f  |   %12.4f  |   %12.4f\n', meanMolecularWeight(mix1), meanMolecularWeight(mix2), meanMolecularWeight(mix3));
    fprintf('(dlV/dlp)T [-] |                 |   %12.4f  |   %12.4f\n', mix2.dVdp_T, mix3.dVdp_T);
    fprintf('(dlV/dlT)p [-] |                 |   %12.4f  |   %12.4f\n', mix2.dVdT_p, mix3.dVdT_p);
    fprintf('cp [kJ/(kg-K)] |   %12.4f  |   %12.4f  |   %12.4f\n', cp_mass(mix1), cp_mass(mix2), cp_mass(mix3));
    fprintf('gamma [-]      |   %12.4f  |   %12.4f  |   %12.4f\n', adiabaticIndex(mix1), adiabaticIndex_sound(mix2), adiabaticIndex_sound(mix3));
    fprintf('sound vel [m/s]|   %12.4f  |   %12.4f  |   %12.4f\n', soundspeed(mix1), soundspeed(mix2), soundspeed(mix3));
    fprintf('u [m/s]        |   %12.4f  |   %12.4f  |   %12.4f\n', velocity_relative(mix1), mix2.v_shock, mix3.v_shock);
    fprintf('Mach number [-]|   %12.4f  |   %12.4f  |   %12.4f\n', velocity_relative(mix1)/soundspeed(mix1), mix2.v_shock/soundspeed(mix2), mix3.v_shock/soundspeed(mix3));
    if contains(ProblemType, '_OBLIQUE') || contains(ProblemType, '_POLAR')
        fprintf('------------------------------------------------------------------------\n');
        fprintf('PARAMETERS\n');
        fprintf('min wave  [deg]|                 |   %12.4f  |   %12.4f\n', mix2.beta_min, mix3.beta_min);
        if contains(ProblemType, '_OBLIQUE')
            fprintf('wave angle[deg]|                 |   %12.4f  |   %12.4f\n', mix2.beta, mix3.beta);
            fprintf('deflection[deg]|                 |   %12.4f  |   %12.4f\n', mix2.theta, mix3.theta);
        else
            fprintf('max def.  [deg]|                 |   %12.4f  |   %12.4f\n', mix2.theta_max);
            fprintf('sonic def.[deg]|                 |   %12.4f  |   %12.4f\n', mix2.theta_sonic);
        end
    elseif strcmpi(ProblemType, 'ROCKET')
        fprintf('------------------------------------------------------------------------\n');
        fprintf('PERFORMANCE PARAMETERS\n');    
        fprintf('CSTAR [m/s]    |                 |   %12.4f  |\n', mix3.cstar);
        fprintf('CF [-]         |                 |   %12.4f  |\n', mix3.cf);
        fprintf('Ivac [m/s]     |                 |   %12.4f  |\n', mix3.I_vac);
        fprintf('Isp  [m/s]     |                 |   %12.4f  |\n', mix3.I_sp);
    end
    fprintf('------------------------------------------------------------------------\n');

    if contains(ProblemType, '_R') || contains(ProblemType, '_OBLIQUE')
        fprintf('STATE 1                 Xi [-]\n');
    else
        fprintf('INLET CHAMBER           Xi [-]\n');
    end
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mix1.Xi(:), ind_sort] = sort(mix1.Xi(:), 'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = mix1.Xi>mintol_display;
    Xminor = sum(mix1.Xi(~j));
    for i=1:length(j)
        if j(i)
              fprintf('%-20s %1.4e\n', ListSpecies{ind_sort(i)}, mix1.Xi(i));
        end
    end
    Nminor = length(mix1.Xi) - sum(j);
    s_space_Nminor = char(32 * ones(1, 4 - numel(num2str(Nminor))));
    fprintf('MINORS[+%d] %s     %12.4e\n\n', Nminor, s_space_Nminor, Xminor);
    fprintf('TOTAL            %14.4e\n',sum(mix1.Xi));
    fprintf('------------------------------------------------------------------------\n');
    if contains(ProblemType, '_R')
        fprintf('STATE 2                 Xi [-]\n');
    elseif contains(ProblemType, '_OBLIQUE')
        fprintf('STATE 2-WEAK SHOCK      Xi [-]\n');
    else
        fprintf('OUTLET CHAMBER          Xi [-]\n');
    end
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mix2.Xi(:), ind_sort] = sort(mix2.Xi(:),'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = mix2.Xi>mintol_display;
    Xminor = sum(mix2.Xi(~j));

    for i=1:length(j)
        if j(i)
              fprintf('%-20s %1.4e\n', ListSpecies{ind_sort(i)}, mix2.Xi(i));
        end
    end
    Nminor = length(mix2.Xi) - sum(j);
    s_space_Nminor = char(32 * ones(1, 4 - numel(num2str(Nminor))));
    fprintf('MINORS[+%d] %s     %12.4e\n\n', Nminor, s_space_Nminor, Xminor);
    fprintf('TOTAL            %14.4e\n', sum(mix2.Xi));
    fprintf('------------------------------------------------------------------------\n');
    if contains(ProblemType, '_R')
        fprintf('STATE 3                 Xi [-]\n');
    elseif contains(ProblemType, '_OBLIQUE')
        fprintf('STATE 2-STRONG SHOCK    Xi [-]\n');
    else
        fprintf('THROAT                  Xi [-]\n');
    end
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mix3.Xi(:), ind_sort] = sort(mix3.Xi(:), 'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    j = mix3.Xi>mintol_display;
    Xminor = sum(mix3.Xi(~j));

    for i=1:length(j)
        if j(i)
              fprintf('%-20s %1.4e\n', ListSpecies{ind_sort(i)}, mix3.Xi(i));
        end
    end
    Nminor = length(mix3.Xi)-sum(j);
    s_space_Nminor = char(32 * ones(1, 4 - numel(num2str(Nminor))));
    fprintf('MINORS[+%d] %s     %12.4e\n\n', Nminor, s_space_Nminor, Xminor);
    fprintf('TOTAL            %14.4e\n',sum(mix3.Xi));
    fprintf('--------------------------------------------------------------------------------------\n');
    fprintf('**************************************************************************************\n\n\n');
elseif nargin == 7
    mix1 = varargin{1}; mix2 = varargin{2}; mix3 = varargin{3}; mix4 = varargin{4};
    fprintf('**************************************************************************************\n');
    fprintf('--------------------------------------------------------------------------------------\n');
    fprintf('Problem type: %s  | phi = %4.4f\n',ProblemType, equivalenceRatio(mix1));
    fprintf('--------------------------------------------------------------------------------------\n');
    if contains(ProblemType, '_R')
        fprintf('               |     STATE 1     |     STATE 2     |     STATE 3     |     STATE 4\n');
    elseif contains(ProblemType, '_OBLIQUE')
        fprintf('               |     STATE 1     |     STATE 2     |     STATE 3-W   |     STATE 3-S\n');
    else
        fprintf('               |  INLET CHAMBER  | OUTLET CHAMBER  |     THROAT      |     EXIT\n');
    end
    
    fprintf('T [K]          |   %12.4f  |   %12.4f  |   %12.4f  |   %12.4f\n', temperature(mix1), temperature(mix2), temperature(mix3), temperature(mix4));
    fprintf('p [bar]        |   %12.4f  |   %12.4f  |   %12.4f  |   %12.4f\n', pressure(mix1), pressure(mix2), pressure(mix3), pressure(mix4));
    fprintf('r [kg/m3]      |   %12.4f  |   %12.4f  |   %12.4f  |   %12.4f\n', density(mix1), density(mix2), density(mix3), density(mix4));
    fprintf('h [kJ/kg]      |   %12.4f  |   %12.4f  |   %12.4f  |   %12.4f\n', enthalpy_mass(mix1), enthalpy_mass(mix2), enthalpy_mass(mix3), enthalpy_mass(mix4));
    fprintf('e [kJ/kg]      |   %12.4f  |   %12.4f  |   %12.4f  |   %12.4f\n', intEnergy_mass(mix1), intEnergy_mass(mix2), intEnergy_mass(mix3), intEnergy_mass(mix4));
    fprintf('g [kJ/kg]      |   %12.4f  |   %12.4f  |   %12.4f  |   %12.4f\n', gibbs_mass(mix1), gibbs_mass(mix2), gibbs_mass(mix3), gibbs_mass(mix4));
    fprintf('s [kJ/(kg-K)]  |   %12.4f  |   %12.4f  |   %12.4f  |   %12.4f\n', entropy_mass(mix1), entropy_mass(mix2), entropy_mass(mix3), entropy_mass(mix4));
    fprintf('W [g/mol]      |   %12.4f  |   %12.4f  |   %12.4f  |   %12.4f\n', meanMolecularWeight(mix1), meanMolecularWeight(mix2), meanMolecularWeight(mix3), meanMolecularWeight(mix4));
    fprintf('(dlV/dlp)T [-] |                 |   %12.4f  |   %12.4f  |   %12.4f\n', mix2.dVdp_T, mix3.dVdp_T, mix4.dVdp_T);
    fprintf('(dlV/dlT)p [-] |                 |   %12.4f  |   %12.4f  |   %12.4f\n', mix2.dVdT_p, mix3.dVdT_p, mix4.dVdT_p);
    fprintf('cp [kJ/(kg-K)] |   %12.4f  |   %12.4f  |   %12.4f  |   %12.4f\n', cp_mass(mix1), cp_mass(mix2), cp_mass(mix3), cp_mass(mix4));
    fprintf('gamma [-]      |   %12.4f  |   %12.4f  |   %12.4f  |   %12.4f\n', adiabaticIndex(mix1), adiabaticIndex_sound(mix2), adiabaticIndex_sound(mix3), adiabaticIndex_sound(mix4));
    fprintf('sound vel [m/s]|   %12.4f  |   %12.4f  |   %12.4f  |   %12.4f\n', soundspeed(mix1), soundspeed(mix2), soundspeed(mix3), soundspeed(mix4));
    fprintf('u [m/s]        |   %12.4f  |   %12.4f  |   %12.4f  |   %12.4f\n', velocity_relative(mix1), mix2.v_shock, mix3.v_shock, mix4.v_shock);
    fprintf('Mach number [-]|   %12.4f  |   %12.4f  |   %12.4f  |   %12.4f\n', velocity_relative(mix1)/soundspeed(mix1), mix2.v_shock/soundspeed(mix2), mix3.v_shock/soundspeed(mix3), mix4.v_shock/soundspeed(mix4));
    if contains(ProblemType, '_OBLIQUE') || contains(ProblemType, '_POLAR')
        fprintf('--------------------------------------------------------------------------------------\n');
        fprintf('PARAMETERS\n');
        fprintf('min wave  [deg]|                 |   %12.4f  |   %12.4f  |   %12.4f\n', mix2.beta_min, mix3.beta_min, mix4.beta_min);
        if contains(ProblemType, '_OBLIQUE')
            fprintf('wave angle[deg]|                 |   %12.4f  |   %12.4f  |   %12.4f\n', mix2.beta, mix3.beta, mix4.beta);
            fprintf('deflection[deg]|                 |   %12.4f  |   %12.4f  |   %12.4f\n', mix2.theta, mix3.theta, mix4.theta);
        else
            fprintf('max def.  [deg]|                 |   %12.4f  |   %12.4f  |   %12.4f\n', mix2.theta_max, mix3.theta_max, mix4.theta_max);
            fprintf('sonic def.[deg]|                 |   %12.4f  |   %12.4f  |   %12.4f\n', mix2.theta_sonic, mix3.theta_sonic, mix4.theta_sonic);
        end
    elseif strcmpi(ProblemType, 'ROCKET')
        fprintf('--------------------------------------------------------------------------------------\n');
        fprintf('PERFORMANCE PARAMETERS\n');    
        fprintf('CSTAR [m/s]    |                 |   %12.4f  |   %12.4f  |\n', mix3.cstar, mix4.cstar);
        fprintf('CF [-]         |                 |   %12.4f  |   %12.4f  |\n', mix3.cf, mix4.cf);
        fprintf('Ivac [m/s]     |                 |   %12.4f  |   %12.4f  |\n', mix3.I_vac, mix4.I_vac);
        fprintf('Isp  [m/s]     |                 |   %12.4f  |   %12.4f  |\n', mix3.I_sp, mix4.I_sp);
    end
    fprintf('--------------------------------------------------------------------------------------\n');

    if contains(ProblemType, '_R') || contains(ProblemType, '_OBLIQUE')
        fprintf('STATE 1                 Xi [-]\n');
    else
        fprintf('INLET CHAMBER           Xi [-]\n');
    end
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mix1.Xi(:), ind_sort] = sort(mix1.Xi(:), 'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = mix1.Xi>mintol_display;
    Xminor = sum(mix1.Xi(~j));
    for i=1:length(j)
        if j(i)
              fprintf('%-20s %1.4e\n', ListSpecies{ind_sort(i)}, mix1.Xi(i));
        end
    end
    Nminor = length(mix1.Xi) - sum(j);
    s_space_Nminor = char(32 * ones(1, 4 - numel(num2str(Nminor))));
    fprintf('MINORS[+%d] %s     %12.4e\n\n', Nminor, s_space_Nminor, Xminor);
    fprintf('TOTAL            %14.4e\n',sum(mix1.Xi));
    fprintf('--------------------------------------------------------------------------------------\n');
    if contains(ProblemType, '_R')
        fprintf('STATE 2                 Xi [-]\n');
    else
        fprintf('OUTLET CHAMBER          Xi [-]\n');
    end
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mix2.Xi(:), ind_sort] = sort(mix2.Xi(:),'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = mix2.Xi>mintol_display;
    Xminor = sum(mix2.Xi(~j));

    for i=1:length(j)
        if j(i)
              fprintf('%-20s %1.4e\n', ListSpecies{ind_sort(i)}, mix2.Xi(i));
        end
    end
    Nminor = length(mix2.Xi) - sum(j);
    s_space_Nminor = char(32 * ones(1, 4 - numel(num2str(Nminor))));
    fprintf('MINORS[+%d] %s     %12.4e\n\n', Nminor, s_space_Nminor, Xminor);
    fprintf('TOTAL            %14.4e\n', sum(mix2.Xi));
    fprintf('--------------------------------------------------------------------------------------\n');
    if contains(ProblemType, '_R')
        fprintf('STATE 3                 Xi [-]\n');
    elseif contains(ProblemType, '_OBLIQUE')
        fprintf('STATE 3-WEAK SHOCK      Xi [-]\n');
    else
        fprintf('THROAT                  Xi [-]\n');
    end
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mix3.Xi(:), ind_sort] = sort(mix3.Xi(:), 'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    j = mix3.Xi>mintol_display;
    Xminor = sum(mix3.Xi(~j));

    for i=1:length(j)
        if j(i)
              fprintf('%-20s %1.4e\n', ListSpecies{ind_sort(i)}, mix3.Xi(i));
        end
    end
    Nminor = length(mix3.Xi)-sum(j);
    s_space_Nminor = char(32 * ones(1, 4 - numel(num2str(Nminor))));
    fprintf('MINORS[+%d] %s     %12.4e\n\n', Nminor, s_space_Nminor, Xminor);
    fprintf('TOTAL            %14.4e\n',sum(mix3.Xi));
    fprintf('--------------------------------------------------------------------------------------\n');
    if contains(ProblemType, '_R')
        fprintf('STATE 4                 Xi [-]\n');
    elseif contains(ProblemType, '_OBLIQUE')
        fprintf('STATE 3-STRONG SHOCK    Xi [-]\n');
    else
        fprintf('EXIT                    Xi [-]\n');
    end
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mix4.Xi(:), ind_sort] = sort(mix4.Xi(:),'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = mix4.Xi>mintol_display;
    Xminor = sum(mix4.Xi(~j));

    for i=1:length(j)
        if j(i)
              fprintf('%-20s %1.4e\n', ListSpecies{ind_sort(i)}, mix4.Xi(i));
        end
    end
    Nminor = length(mix4.Xi) - sum(j);
    s_space_Nminor = char(32 * ones(1, 4 - numel(num2str(Nminor))));
    fprintf('MINORS[+%d] %s     %12.4e\n\n', Nminor, s_space_Nminor, Xminor);
    fprintf('TOTAL            %14.4e\n', sum(mix4.Xi));
    fprintf('--------------------------------------------------------------------------------------\n');
    fprintf('**************************************************************************************\n\n\n');
else
    error('Function displayresults - Not enough arguments')
end