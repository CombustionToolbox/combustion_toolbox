function [alpha, beta, gamma, omega] = set_alpha_beta_gamma_omega(Element_matrix,ind_C,ind_H,ind_O,ind_N)

alpha = 0; if ~isempty(find(Element_matrix(1,:)==ind_C, 1)), alpha = Element_matrix(2,Element_matrix(1,:)==ind_C); end % Number of C-atoms in minor_products{n}
beta  = 0; if ~isempty(find(Element_matrix(1,:)==ind_H, 1)), beta  = Element_matrix(2,Element_matrix(1,:)==ind_H); end % Number of H-atoms in minor_products{n}
gamma = 0; if ~isempty(find(Element_matrix(1,:)==ind_O, 1)), gamma = Element_matrix(2,Element_matrix(1,:)==ind_O); end % Number of N-atoms in minor_products{n}
omega = 0; if ~isempty(find(Element_matrix(1,:)==ind_N, 1)), omega = Element_matrix(2,Element_matrix(1,:)==ind_N); end % Number of O-atoms in minor_products{n}

% alpha = 0; if ~isempty(find(Element_matrix(1,:)==6)), alpha = Element_matrix(2,find(Element_matrix(1,:)==6)); end % Number of C-atoms in minor_products{n}
% beta  = 0; if ~isempty(find(Element_matrix(1,:)==1)), beta  = Element_matrix(2,find(Element_matrix(1,:)==1)); end % Number of H-atoms in minor_products{n}
% gamma = 0; if ~isempty(find(Element_matrix(1,:)==8)), gamma = Element_matrix(2,find(Element_matrix(1,:)==8)); end % Number of N-atoms in minor_products{n}
% omega = 0; if ~isempty(find(Element_matrix(1,:)==7)), omega = Element_matrix(2,find(Element_matrix(1,:)==7)); end % Number of O-atoms in minor_products{n}
