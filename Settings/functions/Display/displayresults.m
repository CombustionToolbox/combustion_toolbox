function displayresults(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   varargin = (strR,str2,strP) or varargin = (strR,strP)
%   mix1  = Prop. state 1 (phi,species,...)
%   mix2  = Prop. state 2 (phi,species,...)
%   mix3  = Prop. state 3 (phi,species,...) - only reflected shocks
% OUTPUT:
%   results on command window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% help displayresults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ProblemType = varargin{end-2};
mintol_display = varargin{end-1};
NameSpecies = varargin{end};

if nargin == 5
    mix1 = varargin{1}; mix2 = varargin{2};
    
    fprintf('***********************************************************\n');
    fprintf('-----------------------------------------------------------\n');
    fprintf('Problem type: %s  | phi = %4.3f\n',ProblemType, equivalenceRatio(mix1));
    fprintf('-----------------------------------------------------------\n');
    if strcmpi(ProblemType,'SHOCK_I') || contains(ProblemType,'DET')
        fprintf('\t\t\t  |   STATE 1   |      STATE 2\n');
    else
        fprintf('\t\t\t  |   REACTANTS |     PRODUCTS\n');
    end
    fprintf('T [K]         |\t %10.3f\t|\t%10.3f\n', temperature(mix1), temperature(mix2));
    fprintf('p [bar]       |\t %10.3f\t|\t%10.3f\n', pressure(mix1), pressure(mix2));
    fprintf('r [kg/m3]     |\t %10.3f\t|\t%10.3f\n', density(mix1), density(mix2));
    fprintf('h [kJ/kg]     |\t %10.3f\t|\t%10.3f\n', enthalpy_mass(mix1), enthalpy_mass(mix2));
    fprintf('e [kJ/kg]     |\t %10.3f\t|\t%10.3f\n', intEnergy_mass(mix1), intEnergy_mass(mix2));
    fprintf('g [kJ/kg]     |\t %10.3f\t|\t%10.3f\n', gibbs_mass(mix1), gibbs_mass(mix2));
    fprintf('s [kJ/(kg-K)] |\t %10.3f\t|\t%10.3f\n', entropy_mass(mix1), entropy_mass(mix2));
    fprintf('cp [kJ/(kg-K)]|\t %10.3f\t|\t%10.3f\n', cp_mass(mix1), cp_mass(mix2));
    fprintf('gamma [-]     |\t %10.3f\t|\t%10.3f\n', adiabaticIndex(mix1), adiabaticIndex(mix2));
    if strcmpi(ProblemType,'SHOCK_I')||strcmpi(ProblemType,'SHOCK_R')||contains(ProblemType,'DET')
        fprintf('u [m/s]       |\t %10.3f\t|\t%10.3f\n', velocity_relative(mix1), velocity_relative(mix2));
    end
    fprintf('-----------------------------------------------------------\n');
    fprintf('REACTANTS\t\t Xi [-]\n');
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mix1.Xi(:),ind_sort] = sort(mix1.Xi(:),'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = mix1.Xi>mintol_display;
    Xminor = sum(mix1.Xi(~j));
    for i=1:length(j)
        if j(i)
            fprintf('%-12s \t%11.4e\n',NameSpecies{ind_sort(i)},mix1.Xi(i));
%             fprintf('%-10s \t%11.4e\n',NameSpecies{i},strR.Xi(i))
        end
    end
    fprintf('MINORS[+%d]   %12.4e\n\n',length(mix1.Xi)-sum(j),Xminor);
    fprintf('TOTAL  \t\t %14.4e\n',sum(mix1.Xi));
    fprintf('-----------------------------------------------------------\n');
    fprintf('PRODUCTS\t\t Xi [-]\n');
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mix2.Xi(:),ind_sort] = sort(mix2.Xi(:),'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = mix2.Xi>mintol_display;
    Xminor = sum(mix2.Xi(~j));
    for i=1:length(j)
        if j(i)
            fprintf('%-12s \t%11.4e\n',NameSpecies{ind_sort(i)},mix2.Xi(i));
%             fprintf('%-10s \t%11.4e\n',NameSpecies{i},strP.Xi(i))
        end
    end
    fprintf('MINORS[+%d]   %12.4e\n\n',length(mix2.Xi)-sum(j),Xminor);
    fprintf('TOTAL  \t\t %14.4e\n',sum(mix2.Xi));
    fprintf('-----------------------------------------------------------\n');
    fprintf('***********************************************************\n\n\n');
elseif nargin == 6
    mix1 = varargin{1}; mix2 = varargin{2}; mix3 = varargin{3};
    fprintf('***********************************************************\n');
    fprintf('-----------------------------------------------------------\n');
    fprintf('Problem type: %s  | phi = %4.3f\n',ProblemType, equivalenceRatio(mix1));
    fprintf('-----------------------------------------------------------\n');
    fprintf('\t\t\t |    STATE 1 \t|    STATE 2 \t|    STATE 3\n');
    fprintf('T [K]        |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n', temperature(mix1), temperature(mix2), temperature(mix3));
    fprintf('p [bar]      |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n', pressure(mix1), pressure(mix2), pressure(mix3));
    fprintf('r [kg/m3]    |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n', density(mix1), density(mix2), density(mix3));
    fprintf('h [kJ/kg]    |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n', enthalpy_mass(mix1), enthalpy_mass(mix2), enthalpy_mass(mix3));
    fprintf('e [kJ/kg]    |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n', intEnergy_mass(mix1), intEnergy_mass(mix2), intEnergy_mass(mix3));
    fprintf('g [kJ/kg]    |\t %10.3f\t|\t%10.3f\n|\t%10.3f\n', gibbs_mass(mix1), gibbs_mass(mix2),  gibbs_mass(mix3));
    fprintf('s [kJ/(kg-K)]|\t %10.3f\t|\t%10.3f\t|\t%10.3f\n', entropy_mass(mix1), entropy_mass(mix2), entropy_mass(mix3));
    fprintf('u [m/s]      |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n', velocity_relative(mix1), velocity_relative(mix2), velocity_relative(mix3));
    fprintf('-----------------------------------------------------------\n');
    fprintf('STATE 1\t\t\t Xi [-]\n');
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mix1.Xi(:),ind_sort] = sort(mix1.Xi(:),'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = mix1.Xi>mintol_display;
    Xminor = sum(mix1.Xi(~j));
    for i=1:length(j)
        if j(i)
              fprintf('%-12s \t%11.4e\n',NameSpecies{ind_sort(i)},mix1.Xi(i));
        end
    end
    fprintf('MINORS[+%d]   %12.4e\n\n',length(mix1.Xi)-sum(j),Xminor);
    fprintf('TOTAL  \t\t %14.4e\n',sum(mix1.Xi));
    fprintf('-----------------------------------------------------------\n');
    fprintf('STATE 2\t\t\t Xi [-]\n');
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mix2.Xi(:),ind_sort] = sort(mix2.Xi(:),'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = mix2.Xi>mintol_display;
    Xminor = sum(mix2.Xi(~j));

    for i=1:length(j)
        if j(i)
%             fprintf('%-10s \t%11.4e\n',NameSpecies{i},str2.Xi(i))
              fprintf('%-12s \t%11.4e\n',NameSpecies{ind_sort(i)},mix2.Xi(i));
        end
    end
    fprintf('MINORS[+%d]   %12.4e\n\n',length(mix2.Xi)-sum(j),Xminor);
    fprintf('TOTAL  \t\t %14.4e\n',sum(mix2.Xi));
    fprintf('-----------------------------------------------------------\n');
    fprintf('STATE 3\t\t\t Xi [-]\n');
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mix3.Xi(:),ind_sort] = sort(mix3.Xi(:),'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    j = mix3.Xi>mintol_display;
    Xminor = sum(mix3.Xi(~j));

    for i=1:length(j)
        if j(i)
%             fprintf('%-10s \t%11.4e\n',NameSpecies{i},str3.Xi(i))
              fprintf('%-12s \t%11.4e\n',NameSpecies{ind_sort(i)},mix3.Xi(i));
        end
    end
    fprintf('MINORS[+%d]   %12.4e\n\n',length(mix3.Xi)-sum(j),Xminor);
    fprintf('TOTAL  \t\t %14.4e\n',sum(mix3.Xi));
    fprintf('-----------------------------------------------------------\n');
    fprintf('***********************************************************\n\n\n');
else
    error('Function displayresults - Not enough arguments')
end