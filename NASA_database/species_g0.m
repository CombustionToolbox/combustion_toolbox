function g0i = species_g0(Species,T,strThProp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculates the gibbs free energy (in kJ/mol) at the specified temperature 
% T for the chemical species
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    g0i = strThProp.(Species).g0curve(T)/1000;
catch
    g0i = strThProp.(Species).g0/1000;
end
%  tic, for T=300:100:3000 ; species_g0_new(strThProp,'CO',1800), end, toc

% Tref = 273.15;
% 
% Element_matrix = set_element_matrix(strThProp.(Species).txFormula,strThProp.('CO2').elements); % sets Element_matrix matrix
% set_reference_form_of_elements_with_T_intervals; % sets Reference_form_of_elements_with_T_intervals list
% 
% 
% Species(strfind(Species,'Al')+1)='L';
% Species(strfind(Species,'Cl')+1)='L';
% Species(strfind(Species,'Tl')+1)='L';
% Species(strfind(Species,'Fl')+1)='L';
% Species_with_parenthesis = Species;
% n_open_parenthesis = detect_location_of_phase_specifier(Species_with_parenthesis);
% 
% for i = 1:strThProp.(Species).ctTInt
%     if (T >= strThProp.(Species).tRange{i}(1)) && (T <= strThProp.(Species).tRange{i}(2))
%         tInterval = i;
%     end
% end
% if ~exist('tInterval','var')
%     g0i = strThProp.(Species).g0curve(T)/1000; 
%     return
% end
% R0 = 8.3144598;
% Cp0 = R0 *      sum(strThProp.(Species).a{tInterval} .* T.^strThProp.(Species).tExponents{tInterval});
% H0  = R0 * T * (sum(strThProp.(Species).a{tInterval} .* T.^strThProp.(Species).tExponents{tInterval} .* [-1   log(T) 1      1/2 1/3 1/4 1/5 0]) + strThProp.(Species).b{tInterval}(1)/T);
% S0  = R0 *     (sum(strThProp.(Species).a{tInterval} .* T.^strThProp.(Species).tExponents{tInterval} .* [-1/2 -1     log(T) 1   1/2 1/3 1/4 0]) + strThProp.(Species).b{tInterval}(2)  );
%     
% [iRE, REname] = isRefElm(Reference_form_of_elements_with_T_intervals,Species(1:n_open_parenthesis-1),T);
% if (~iRE)
%     GP = H0 - T.*S0;
%     
%     GR = zeros(1,size(Element_matrix,2));
%     for i = 1:size(Element_matrix,2)
%         nu_i = Element_matrix(2,i);
%         [~, REname_i] = isRefElm(Reference_form_of_elements_with_T_intervals,strThProp.('CO2').elements{Element_matrix(1,i)},T);
% %             [~, REname_i] = isRefElm(Reference_form_of_elements_with_T_intervals,upper(Elements{Element_matrix(1,i)}),T);
% %             [~, ~, ~, ~, ~, H0_i, ~, ~, S0_i, ~] = SpeciesThermProp(strMaster,REname_i,T,'molar',0);
%         [~, ~, ~, ~, ~, H0_i, ~, ~, S0_i, ~] = SpeciesThermProp(strThProp,REname_i,T,'molar',0);
%         GR(i) = nu_i*(H0_i - T.*S0_i);
%         if any(Element_matrix(1,i)==[1, 7, 8, 9, 17, 35]), GR(i) = GR(i)/2; end
%     end
%     GR = sum(GR);
%     g0i = (GP-GR)/1000;
% else
%     g0i = 0;
% end