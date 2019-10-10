function Species_matrix = set_species_matrix(strMaster,varargin)
global MassOrMolar
if nargin < 2
    disp(['function set_species - Not enough arguments.'])
    return
end

T = 298.15;

Species_matrix = [];
if mod(length(varargin),2)==0
    for n = 1:length(varargin)/2
        Species = varargin{2*n-1};
        if strcmp(Species,'T')
            T = varargin{2*n};
        end
    end    
    for n = 1:length(varargin)/2
        Species = varargin{2*n-1};
        nu      = varargin{2*n};
        if ~strcmp(Species,'T')
            [txFormula, ~, ~, ~, ~, ~, ~, ~, ~, ~] = SpeciesThermProp(strMaster,Species,T,MassOrMolar);
            Element_matrix = set_element_matrix(txFormula);
            Element_matrix(2,:) = nu*Element_matrix(2,:);
            if n == 1
                Species_matrix = Element_matrix;
            else
                for i = 1:size(Element_matrix,2)
                    j = find(Element_matrix(1,i)==Species_matrix(1,:));
                    if isempty(j)
                        AUX_Element_matrix = [Element_matrix(1,i);zeros(n-1,1);Element_matrix(2,i)];
                        Species_matrix = [[Species_matrix; zeros(1,size(Species_matrix,2))], AUX_Element_matrix];
                    else
                        Species_matrix(n+1,j) = Element_matrix(2,i);
                    end
                end
            end    
        end
    end
    Species_matrix = sortrows(Species_matrix',1)';
    n_elements = size(Species_matrix,2);
    for n = 1:length(varargin)/2
        Species = varargin{2*n-1};
        nu      = varargin{2*n};
        if ~strcmp(Species,'T')
%             [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = SpeciesThermProp(strMaster,Species,T,MassOrMolar);
            [~, ~, Cp0, Cv0, ~, H0, ~, E0, ~, ~] = SpeciesThermProp(strMaster,Species,T,MassOrMolar);
            Species_matrix(n+1,n_elements+1:4) = nu*[Cp0, Cv0, H0/1000, E0/1000];
        end
    end
end

Species_matrix = [Species_matrix; sum(Species_matrix(2:end,:))];
