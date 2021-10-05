function displayresults(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY RESULTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   varargin = (strR,str2,strP) or varargin = (strR,strP)
%   strR  = Prop. of reactives (phi,species,...)
%   strP  = Prop. of products (phi,species,...)
%   str2  = Prop. of state 2 (phi,species,...)
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
    fprintf('Problem type: %s  | phi = %4.3f\n',ProblemType,mix1.phi);
    fprintf('-----------------------------------------------------------\n');
    if strcmpi(ProblemType,'SHOCK_I') || contains(ProblemType,'DET')
        fprintf('\t\t\t  |   STATE 1   |      STATE 2\n');
    else
        fprintf('\t\t\t  |   REACTANTS |     PRODUCTS\n');
    end
    fprintf('T [K]         |\t %10.3f\t|\t%10.3f\n',mix1.T,mix2.T);
    fprintf('p [bar]       |\t %10.3f\t|\t%10.3f\n',mix1.p,mix2.p);
    fprintf('r [kg/m3]     |\t %10.3f\t|\t%10.3f\n',mix1.rho,mix2.rho);
    fprintf('h [kJ/kg]     |\t %10.3f\t|\t%10.3f\n',mix1.h/mix1.mi,mix2.h/mix2.mi);
    fprintf('e [kJ/kg]     |\t %10.3f\t|\t%10.3f\n',mix1.e/mix1.mi,mix2.e/mix2.mi);
    fprintf('s [kJ/(kg-K)] |\t %10.3f\t|\t%10.3f\n',mix1.S,mix2.S);
    fprintf('cp [kJ/(kg-K)]|\t %10.3f\t|\t%10.3f\n',mix1.cP/mix1.mi*1e-3,mix2.cP/mix2.mi*1e-3);
    fprintf('gamma [-]     |\t %10.3f\t|\t%10.3f\n',mix1.cP/mix1.cV,mix2.cP/mix2.cV);
    if strcmpi(ProblemType,'SHOCK_I')||strcmpi(ProblemType,'SHOCK_R')||contains(ProblemType,'DET')
        fprintf('u [m/s]       |\t %10.3f\t|\t%10.3f\n',mix1.u,mix2.u);
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
    fprintf('Problem type: %s  | phi = %4.3f\n',ProblemType,mix1.phi);
    fprintf('-----------------------------------------------------------\n');
    fprintf('\t\t\t |    STATE 1 \t|    STATE 2 \t|    STATE 3\n');
    fprintf('T [K]        |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n',mix1.T,mix2.T,mix3.T);
    fprintf('p [bar]      |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n',mix1.p,mix2.p,mix3.p);
    fprintf('r [kg/m3]    |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n',mix1.rho,mix2.rho,mix3.rho);
    fprintf('h [kJ/kg]    |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n',mix1.h/mix1.mi,mix2.h/mix2.mi,mix3.h/mix3.mi);
    fprintf('e [kJ/kg]    |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n',mix1.e/mix1.mi,mix2.e/mix2.mi,mix3.e/mix3.mi);
    fprintf('s [kJ/(kg-K)]|\t %10.3f\t|\t%10.3f\t|\t%10.3f\n',mix1.S,mix2.S,mix3.S);
    fprintf('u [m/s]      |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n',mix1.u,mix2.u,mix3.u);
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