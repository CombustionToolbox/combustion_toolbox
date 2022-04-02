function element_matrix = set_element_matrix(txFormula, elements)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE element MATRIX OF SPECIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   txFormula = 
% OUTPUT:
%   element_matrix = element matrix of species
% 
% EXAMPLE: for CO2
%
% element_matrix =
% 
%      7     9
%      1     2
%
% That is, the species contains 1 atom of element 7 (C) and 2 atoms of
% element 9 (O)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% help set_element_matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

element_1 = txFormula(1:2);
if strcmp(element_1, '  ')
else
    if strcmp(element_1(2), ' ')
        element_1 = element_1(1);
    end
    element_matrix(1,1) = find(strcmpi(elements, element_1));
    element_matrix(2,1) = sscanf(txFormula(3:8), '%f');
end

element_2 = txFormula(9:10);
if strcmp(element_2, '  ')
else
    if strcmp(element_2(2), ' ')
        element_2 = element_2(1); 
    end
    element_matrix(1,2) = find(strcmpi(elements, element_2));
    element_matrix(2,2) = sscanf(txFormula(11:16), '%f');
end

element_3 = txFormula(17:18);
if strcmp(element_3, '  ')
else
    if strcmp(element_3(2), ' ')
        element_3 = element_3(1); 
    end
    element_matrix(1,3) = find(strcmpi(elements, element_3));
    element_matrix(2,3) = sscanf(txFormula(19:24), '%f');
end

element_4 = txFormula(25:26);
if strcmp(element_4, '  ')
else
    if strcmp(element_4(2), ' ')
        element_4 = element_4(1); 
    end
    element_matrix(1,4) = find(strcmpi(elements, element_4));
    element_matrix(2,4) = sscanf(txFormula(27:32), '%f');
end

element_5 = txFormula(33:34);
if strcmp(element_5, '  ')
else
    if strcmp(element_5(2), ' ')
        element_5 = element_5(1); 
    end
    element_matrix(1,5) = find(strcmpi(elements, element_5));
    element_matrix(2,5) = sscanf(txFormula(35:40), '%f');
end
