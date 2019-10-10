function displayresults_m2i(varargin)
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
    fprintf('REACTANTS\t\t Xi [-]\n');
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [strR.Xi(:),ind_sort] = sort(strR.Xi(:)*strR.N,'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = strR.Xi>mintol_display;
    Xminor = sum(strR.Xi(~j));
    for i=1:length(j)
        if j(i)
            fprintf('%-12s \t%11.4f\n',NameSpecies{ind_sort(i)},strR.Xi(i));
%             fprintf('%-10s \t%11.4f\n',NameSpecies{i},strR.Xi(i))
        end
    end
    fprintf('MINORS[+%d]   %12.4f\n\n',length(strR.Xi)-sum(j),Xminor);
    fprintf('TOTAL  \t\t %14.4f\n',sum(strR.Xi));
    fprintf('-----------------------------------------------------------\n');
    fprintf('PRODUCTS\t\t Xi [-]\n');
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [strP.Xi(:),ind_sort] = sort(strP.Xi(:)*strP.N,'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = strP.Xi>mintol_display;
    Xminor = sum(strP.Xi(~j));
    for i=1:length(j)
        if j(i)
            fprintf('%-12s \t%11.4f\n',NameSpecies{ind_sort(i)},strP.Xi(i));
%             fprintf('%-10s \t%11.4f\n',NameSpecies{i},strP.Xi(i))
        end
    end
    fprintf('MINORS[+%d]   %12.4f\n\n',length(strP.Xi)-sum(j),Xminor);
    fprintf('TOTAL  \t\t %14.4f\n',sum(strP.Xi));
    fprintf('-----------------------------------------------------------\n');
    fprintf('***********************************************************\n\n\n');
end