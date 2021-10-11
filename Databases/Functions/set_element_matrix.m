function Element_matrix = set_element_matrix(txFormula,Elements)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE ELEMENT MATRIX OF SPECIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   txFormula = 
% OUTPUT:
%   Element_matrix = element matrix of species
% 
% EXAMPLE: for CO2
%
% Element_matrix =
% 
%      6     8
%      1     2
%
% That is, the species contains 1 atom of element 6 (C) and 2 atoms of
% element 8 (O)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% help set_element_matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Element_matrix = zeros(2,5);

Element_1 = txFormula(1:2);
if strcmp(Element_1,'  ')
else
    if strcmp(Element_1(2),' ')
        Element_1 = Element_1(1);
    end
    Element_matrix(1,1) = find(strcmpi(Elements,Element_1));
    Element_matrix(2,1) = sscanf(txFormula(3:8),'%f');
end

Element_2 = txFormula(9:10);
if strcmp(Element_2,'  ')
else
    if strcmp(Element_2(2),' ')
        Element_2 = Element_2(1); 
    end
    Element_matrix(1,2) = find(strcmpi(Elements,Element_2));
    Element_matrix(2,2) = sscanf(txFormula(11:16),'%f');
end

Element_3 = txFormula(17:18);
if strcmp(Element_3,'  ')
else
    if strcmp(Element_3(2),' ')
        Element_3 = Element_3(1); 
    end
    Element_matrix(1,3) = find(strcmpi(Elements,Element_3));
    Element_matrix(2,3) = sscanf(txFormula(19:24),'%f');
end

Element_4 = txFormula(25:26);
if strcmp(Element_4,'  ')
else
    if strcmp(Element_4(2),' ')
        Element_4 = Element_4(1); 
    end
    Element_matrix(1,4) = find(strcmpi(Elements,Element_4));
    Element_matrix(2,4) = sscanf(txFormula(27:32),'%f');
end

Element_5 = txFormula(33:34);
if strcmp(Element_5,'  ')
else
    if strcmp(Element_5(2),' ')
        Element_5 = Element_5(1); 
    end
    Element_matrix(1,5) = find(strcmpi(Elements,Element_5));
    Element_matrix(2,5) = sscanf(txFormula(35:40),'%f');
end

% Element_matrix = sortrows(Element_matrix',1)';
