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
    strR = varargin{1}; strP = varargin{2};
    fprintf('-----------------------------------------------------------\n');
    fprintf('Problem type: %s  | phi = %4.3f\n',ProblemType,strR.phi);
    fprintf('-----------------------------------------------------------\n');
    if strcmpi(ProblemType,'SHOCK_I') || contains(ProblemType,'DET')
        fprintf('\t\t\t  |   STATE 1   |      STATE 2\n');
    else
        fprintf('\t\t\t  |   REACTANTS |     PRODUCTS\n');
    end
    fprintf('T [K]         |\t %10.3f\t|\t%10.3f\n',strR.T,strP.T);
    fprintf('p [bar]       |\t %10.3f\t|\t%10.3f\n',strR.p,strP.p);
    fprintf('r [kg/m3]     |\t %10.3f\t|\t%10.3f\n',strR.rho,strP.rho);
    fprintf('h [kJ/kg]     |\t %10.3f\t|\t%10.3f\n',strR.h/strR.mi,strP.h/strP.mi);
    fprintf('e [kJ/kg]     |\t %10.3f\t|\t%10.3f\n',strR.e/strR.mi,strP.e/strP.mi);
    fprintf('s [kJ/(kg-K)] |\t %10.3f\t|\t%10.3f\n',strR.S,strP.S);
    fprintf('cp [kJ/(kg-K)]|\t %10.3f\t|\t%10.3f\n',strR.cP/strR.mi*1e-3,strP.cP/strP.mi*1e-3);
    fprintf('gamma [-]     |\t %10.3f\t|\t%10.3f\n',strR.cP/strR.cV,strP.cP/strP.cV);
    if strcmpi(ProblemType,'SHOCK_I')||strcmpi(ProblemType,'SHOCK_R')||contains(ProblemType,'DET')
        fprintf('u [m/s]       |\t %10.3f\t|\t%10.3f\n',strR.u,strP.u);
    end
    fprintf('-----------------------------------------------------------\n');
    fprintf('REACTANTS\t\t Xi [-]\n');
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [strR.Xi(:),ind_sort] = sort(strR.Xi(:),'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = strR.Xi>mintol_display;
    Xminor = sum(strR.Xi(~j));
    for i=1:length(j)
        if j(i)
            fprintf('%-12s \t%11.4e\n',NameSpecies{ind_sort(i)},strR.Xi(i));
%             fprintf('%-10s \t%11.4e\n',NameSpecies{i},strR.Xi(i))
        end
    end
    fprintf('MINORS[+%d]   %12.4e\n\n',length(strR.Xi)-sum(j),Xminor);
    fprintf('TOTAL  \t\t %14.4e\n',sum(strR.Xi));
    fprintf('-----------------------------------------------------------\n');
    fprintf('PRODUCTS\t\t Xi [-]\n');
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [strP.Xi(:),ind_sort] = sort(strP.Xi(:),'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = strP.Xi>mintol_display;
    Xminor = sum(strP.Xi(~j));
    for i=1:length(j)
        if j(i)
            fprintf('%-12s \t%11.4e\n',NameSpecies{ind_sort(i)},strP.Xi(i));
%             fprintf('%-10s \t%11.4e\n',NameSpecies{i},strP.Xi(i))
        end
    end
    fprintf('MINORS[+%d]   %12.4e\n\n',length(strP.Xi)-sum(j),Xminor);
    fprintf('TOTAL  \t\t %14.4e\n',sum(strP.Xi));
    fprintf('-----------------------------------------------------------\n');
    fprintf('***********************************************************\n\n\n');
elseif nargin == 6
    str1 = varargin{1}; str2 = varargin{2}; str3 = varargin{3};
    fprintf('-----------------------------------------------------------\n');
    fprintf('Problem type: %s  | phi = %4.3f\n',ProblemType,str1.phi);
    fprintf('-----------------------------------------------------------\n');
    fprintf('\t\t\t |    STATE 1 \t|    STATE 2 \t|    STATE 3\n');
    fprintf('T [K]        |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n',str1.T,str2.T,str3.T);
    fprintf('p [bar]      |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n',str1.p,str2.p,str3.p);
    fprintf('r [kg/m3]    |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n',str1.rho,str2.rho,str3.rho);
    fprintf('h [kJ/kg]    |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n',str1.h/str1.mi,str2.h/str2.mi,str3.h/str3.mi);
    fprintf('e [kJ/kg]    |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n',str1.e/str1.mi,str2.e/str2.mi,str3.e/str3.mi);
    fprintf('s [kJ/(kg-K)]|\t %10.3f\t|\t%10.3f\t|\t%10.3f\n',str1.S,str2.S,str3.S);
    fprintf('u [m/s]      |\t %10.3f\t|\t%10.3f\t|\t%10.3f\n',str1.u,str2.u,str3.u);
    fprintf('-----------------------------------------------------------\n');
    fprintf('STATE 1\t\t\t Xi [-]\n');
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [str1.Xi(:),ind_sort] = sort(str1.Xi(:),'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = str1.Xi>mintol_display;
    Xminor = sum(str1.Xi(~j));
    for i=1:length(j)
        if j(i)
              fprintf('%-12s \t%11.4e\n',NameSpecies{ind_sort(i)},str1.Xi(i));
        end
    end
    fprintf('MINORS[+%d]   %12.4e\n\n',length(str1.Xi)-sum(j),Xminor);
    fprintf('TOTAL  \t\t %14.4e\n',sum(str1.Xi));
    fprintf('-----------------------------------------------------------\n');
    fprintf('STATE 2\t\t\t Xi [-]\n');
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [str2.Xi(:),ind_sort] = sort(str2.Xi(:),'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = str2.Xi>mintol_display;
    Xminor = sum(str2.Xi(~j));

    for i=1:length(j)
        if j(i)
%             fprintf('%-10s \t%11.4e\n',NameSpecies{i},str2.Xi(i))
              fprintf('%-12s \t%11.4e\n',NameSpecies{ind_sort(i)},str2.Xi(i));
        end
    end
    fprintf('MINORS[+%d]   %12.4e\n\n',length(str2.Xi)-sum(j),Xminor);
    fprintf('TOTAL  \t\t %14.4e\n',sum(str2.Xi));
    fprintf('-----------------------------------------------------------\n');
    fprintf('STATE 3\t\t\t Xi [-]\n');
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [str3.Xi(:),ind_sort] = sort(str3.Xi(:),'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    j = str3.Xi>mintol_display;
    Xminor = sum(str3.Xi(~j));

    for i=1:length(j)
        if j(i)
%             fprintf('%-10s \t%11.4e\n',NameSpecies{i},str3.Xi(i))
              fprintf('%-12s \t%11.4e\n',NameSpecies{ind_sort(i)},str3.Xi(i));
        end
    end
    fprintf('MINORS[+%d]   %12.4e\n\n',length(str3.Xi)-sum(j),Xminor);
    fprintf('TOTAL  \t\t %14.4e\n',sum(str3.Xi));
    fprintf('-----------------------------------------------------------\n');
    fprintf('***********************************************************\n\n\n');
else
    error('Function displayresults - Not enough arguments')
end