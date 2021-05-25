function [P_IC,DeltaNP] = CalculateProductsIC(app, N_CC, phi, pP, TP, vP, phi_c, FLAG_SOOT)

% Abbreviations ---------------------
M0 = app.C.M0.value;
A0 = app.C.A0.value;
E = app.E;
S = app.S;
C = app.C;
M = app.M;
PD = app.PD;
TN = app.TN;
strThProp = app.strThProp;
minors_products = app.M.minors_products;
% -----------------------------------


R0TP = C.R0*TP;
it = 0; itMax = 500; t=true;
% Relaxation/iteration parameters
relax = 0.00007385775+(0.9854897-0.00007385775)/(1+(TP/4058911)^1.817875)^658457.8;
% Number of moles of the major species in the product mixture under the
% assumption of complete combustion (CC), denoted by subscript _0
NCO2_0 = N_CC(S.ind_CO2,1);
NCO_0  = N_CC(S.ind_CO,1);
NH2O_0 = N_CC(S.ind_H2O,1);
NH2_0  = N_CC(S.ind_H2,1);
NO2_0  = N_CC(S.ind_O2,1);
NN2_0  = N_CC(S.ind_N2,1);
NCgr_0 = N_CC(S.ind_Cgr,1);
% Number of C, H, O, N, He, Ar-atoms in the product species
NatomE = sum(N_CC(:,1).*A0);

x = NatomE(E.ind_C); %  6 = find(strcmp(Elements,'C'))
y = NatomE(E.ind_H); %  1 = find(strcmp(Elements,'H'))
z = NatomE(E.ind_O); %  8 = find(strcmp(Elements,'O'))
w = NatomE(E.ind_N); %  7 = find(strcmp(Elements,'N'))

% if (phi > 1) && (x*y ~= 0), relax = relax/2; end

% Initial guess for the number of moles of the major species in the
% product species under the assumption of incomplete combustion (IC)
NCO2   = NCO2_0;
NCO    = NCO_0;
NH2O   = NH2O_0;
NH2    = NH2_0;
NO2    = NO2_0;
NN2    = NN2_0;
NCgr   = NCgr_0;
NHe    = NatomE(E.ind_He); %  2 = find(strcmp(Elements,'He'))
NAr    = NatomE(E.ind_Ar); % 18 = find(strcmp(Elements,'Ar'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial guess for the overall number of moles of gaseous species in the
% product mixture

NP = sum(N_CC(:,1).*(1-N_CC(:,2))); % Sum of num of moles of gases-(1-swt), with swt == condensed phase.
% NP = sum(N_CC(:,1));

if strfind(PD.ProblemType,'P') == 2 % PD.ProblemType = 'TP', 'HP', or 'SP'
    zeta = NP/pP;
elseif strfind(PD.ProblemType,'V') == 2 % PD.ProblemType = 'TV', 'EV', or 'SV'
    zeta = (vP/1000)*1e5/R0TP; % vP is given in l
end

P_IC = M0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATION OF GIBBS FREE ENERGY, CONSTANTS EQUILIBRIUM AND OTHERS
g_CO2 = species_g0_new('CO2',TP,strThProp);
g_CO  = species_g0_new('CO',TP,strThProp);
g_H2O = species_g0_new('H2O',TP,strThProp);
DG0_I    = (g_CO-g_CO2)*1000;
DG0_II   = -g_H2O*1000;
k1  = exp(-DG0_I/R0TP);
k2  = exp(-DG0_II/R0TP);
if M.major_CH4
    g_CH4 = species_g0_new('CH4',TP,strThProp);
    g_C2H2 = species_g0_new('C2H2_acetylene',TP,strThProp);
    g_C6H6 = species_g0_new('C6H6',TP,strThProp);
    % DG0_VIII   = g_CH4*1000;
    
    DG0_VIII   = (species_g0_new('CH3',TP,strThProp)-species_g0_new('H',TP,strThProp)-g_CH4)*1000;
    
    k8  = exp(-DG0_VIII/R0TP);
    
    
    DG0_IX = -(g_C6H6-3*g_C2H2)*1000;
    k9 =  exp(-DG0_IX/R0TP);
    
    DG0_X = -(species_g0_new('H',TP,strThProp)+species_g0_new('CH',TP,strThProp))*1000;
    k10 = exp(-DG0_X/R0TP);
end
if M.major_OH
    g_OH = species_g0_new('OH',TP,strThProp);
    DG0_XI   = g_OH*1000;
    
    k11 = exp(-DG0_XI/R0TP);
end
% if M.minor_C
%     DG0_X = 2*g_CO*1000;
%     k10 = exp(-DG0_X/R0TP);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correction to the initial guess of oxygene moles and overall number of
% moles
NO2    = NO2_0+((NH2O*k2+NCO2*k1)/2)^(2/3)*zeta^(1/3);
NP     = NP+(NO2-NO2_0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if phi<=1 && M.Lminors>0 % case of lean-to-stoichiometric mixtures
    DNfactor_III = 1-(C.beta+2*(C.gamma+C.omega))/4;
    for n = M.Lminors:-1:1
        DG0_III(n) = (species_g0_new(minors_products{n},TP,strThProp)-C.alpha(n)*g_CO2 ...
            -(C.beta(n)/2)*g_H2O)*1000;
    end
    k3 = exp(-DG0_III/R0TP);
elseif phi>1  % case rich mixtures
    if (x == 0) && (y ~= 0) && M.Lminors>0 % if there are only hydrogens (H)
        for n = M.Lminors:-1:1
            DG0_VI(n) = (species_g0_new(minors_products{n},TP,strThProp) ...
                -C.alpha(n) * g_CO2 ...
                -(C.gamma(n)-2*C.alpha(n)) * g_H2O)*1000;
        end
        k6 = exp(-DG0_VI/R0TP);
        DNfactor_VI = 1-C.alpha-C.beta/2-C.omega/2;
    elseif ((x ~= 0) && (y == 0) && M.Lminors>0 && phi < phi_c) && ~FLAG_SOOT% if there are only carbons (C)
        for n = M.Lminors:-1:1
            DG0_V(n) = (species_g0_new(minors_products{n},TP,strThProp) ...
                -(C.gamma(n)-C.alpha(n)-C.beta(n)/2) * g_CO2 ...
                -(C.beta(n)/2) * g_H2O ...
                -(2*C.alpha(n)-C.gamma(n)+C.beta(n)/2) * g_CO)*1000;
        end
        k5 = exp(-DG0_V/R0TP);
        DNfactor_V = 1-C.alpha-C.beta/2-C.omega/2;
    elseif phi < phi_c*TN.factor_c && ~FLAG_SOOT% general case of rich mixtures with hydrogens (H) and carbons (C)
        for n = M.Lminors:-1:1
            DG0_V(n) = (species_g0_new(minors_products{n},TP,strThProp) ...
                -(C.gamma(n)-C.alpha(n)-C.beta(n)/2) * g_CO2 ...
                -(C.beta(n)/2) * g_H2O ...
                -(2*C.alpha(n)-C.gamma(n)+C.beta(n)/2) * g_CO)*1000;
        end
        if M.Lminors>0
            k5 = exp(-DG0_V/R0TP);
            DNfactor_V = 1-C.alpha-C.beta/2-C.omega/2;
        end
        DG0_IV = (g_CO+g_H2O-g_CO2)*1000;
        k4 = exp(-DG0_IV/R0TP);
    elseif phi >= phi_c*TN.factor_c || FLAG_SOOT
        for n = M.Lminors:-1:1
            DG0_V(n) = (species_g0_new(minors_products{n},TP,strThProp) ...
                -(C.gamma(n)-C.alpha(n)-C.beta(n)/2) * g_CO2 ...
                -(C.beta(n)/2) * g_H2O ...
                -(2*C.alpha(n)-C.gamma(n)+C.beta(n)/2) * g_CO)*1000;
        end
        if M.Lminors>0
            k5 = exp(-DG0_V/R0TP);
            DNfactor_V = 1-C.alpha-C.beta/2-C.omega/2;
        end
        
        DG0_IV = (g_CO+g_H2O-g_CO2)*1000;
        k4 = exp(-DG0_IV/R0TP);
        DG0_VII = (g_CO2-2*g_CO)*1000;
        k7 = exp(-DG0_VII/(C.R0*TP));
        mu = k7/zeta;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEMICAL EQUILIBRIUM CALCULATIONS
if phi <= 1 % case of lean-to-stoichiometric mixtures
    [P_IC,DeltaNP] = incomplete_phi_1(P_IC,E,S,M,C);  
elseif phi > 1 && (x == 0) && (y ~= 0) %  case of rich mixtures if there are only hydrogens (H)
    [P_IC,DeltaNP] = incomplete_phi_2(P_IC,E,S,M,C);   
elseif ((x ~= 0) && (y == 0) && phi < phi_c) && ~FLAG_SOOT % if there are only carbons (C)
    [P_IC,DeltaNP] = incomplete_phi_3(P_IC,E,S,M,C);     
elseif phi < phi_c*TN.factor_c && ~FLAG_SOOT % general case of rich mixtures with hydrogens (H) and carbons (C) and without soot
    [P_IC,DeltaNP] = incomplete_phi_4(P_IC,E,S,M,C); 
elseif phi >= phi_c*TN.factor_c || FLAG_SOOT% rich mixtures with soot
    if (x ~= 0) && (y == 0) % if there are only carbons (C)
        [P_IC,DeltaNP] = incomplete_phi_5(P_IC,E,S,M,C);
    else % general case
        [P_IC,DeltaNP] = incomplete_phi_6(P_IC,E,S,M,C);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if it>=itMax && strcmp(PD.ProblemType,"TP")
    fprintf('****************************\n');
	fprintf('** Solution not converged **\n');
    fprintf('** phi   =  %4.2f         **\n',phi);
    fprintf('** Temp  =  %4.2f         **\n',TP);
    fprintf('** Error =  %4.2f%%         **\n',abs(DeltaNP/NP)*100);
    fprintf('** It    =  %4.d          **\n',it);
    fprintf('****************************\n');
end
% fprintf('** It    =  %4.d          **\n',it);
% fprintf('** NCgr    =  %14.4e          **\n',NCgr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NESTED FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHEMICAL EQUILIBRIUM
    function [P_IC,DeltaNP] = incomplete_phi_1(P_IC,E,S,M,C)
        DeltaNP = 1;
        while abs(DeltaNP/NP) > C.tolN && it<itMax
            it = it+1;
            % Initialization of the product matrix for incomplete combustion (IC)
            NP_old   = NP;
            P_IC_old = P_IC;
            P_IC = M0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % In lean combustion the product mixture always contains O2, CO2,
            % H2O, and N2 (in fuel-air combustion). The number of moles of
            % CO and H2 can be calculated from the equilibrium condition for
            % the reactions
            %
            %                  CO2 <-I -> CO+(1/2) O2             [k1]
            %                  H2O <-II-> H2+(1/2) O2             [k2]
            %
            % For the remaining minor species, we calculate the number of
            % moles from the equilibrium condition for the generalized
            % reaction
            %
            % C.alpha CO2+(C.beta/2) H2O+(C.gamma/2-C.alpha-C.beta/4) O2
            %+(C.omega/2) N2 <-III-> C_alpha H_beta O_gamma N_omega [k3]
            %
            % Note that reactions*i*and II are particular forms of the more
            % general reaction III
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            NCO    = NCO2/NO2^(1/2)*zeta^(1/2)*k1;
            NCO_old = P_IC_old(S.ind_CO,1);
            
            if NCO_old ~=0
                NCO = NCO_old+relax*(NCO-NCO_old);
            end
            P_IC(S.ind_CO,1) = NCO;
            
            NH2     = NH2O/NO2^(1/2)*zeta^(1/2)*k2;
            NH2_old = P_IC_old(S.ind_H2,1);
            if NH2_old ~=0
                NH2 = NH2_old+relax*(NH2-NH2_old);
            end
            P_IC(S.ind_H2,1) = NH2;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Determination of the number of moles of the minor species
            % from the equilibrium condition for the above reaction
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if M.Lminors>0
                Ni     = k3.* NCO2.^C.alpha.* NH2O.^(C.beta/2).* NN2.^(C.omega/2) ...
                    .* NO2.^(C.gamma/2-C.alpha-C.beta/4).* zeta.^DNfactor_III;
                %           Ni = exp(log(k3)+ log(NCO2).*C.alpha+ log(NH2O).*(C.beta/2)+ log(NN2+1).*(C.omega/2) ...
                %                +log(NO2).*(C.gamma/2-C.alpha-C.beta/4)+ log(zeta).*DNfactor_III);
                %             for n =M.Lminors:-1:1
                %                 Ni_old = P_IC_old(M.ind_minor(n),1);
                %                 if Ni_old ~=0
                %                     Ni(n) = Ni_old+relax*(Ni(n)-Ni_old);
                %                 end
                %                 P_IC(M.ind_minor(n),[strThProp.(minors_products{n}).Element_matrix(1,:),1]) = [Ni(n)*strThProp.(minors_products{n}).Element_matrix(2,:),Ni(n)];
                %             end
%                 if M.major_OH && DeltaNP
%                     Ni(M.ind_m_OH) = sqrt(NH2*NO2*k11*zeta^(-3/2));
%                 end
                Ni(Ni>NP_old) = 0.75*Ni(Ni>NP_old);
                Ni_old = P_IC_old(M.ind_minor,1)';
                aux = find(Ni_old~=0);
                if ~isempty(aux)
                    Ni(aux) = Ni_old(aux)+relax*(Ni(aux)-Ni_old(aux));
                end
                P_IC(M.ind_minor,1) = Ni;
                % disp(sortrows([cell2table(S.NameSpecies(P_IC(:,1)>0.000001), 'VariableNames', {'Species'}), array2table(P_IC(P_IC(:,1)>0.000001), 'VariableNames', {'N'})], 2, 'descend'))
            end
        % disp([sum(P_IC(:,1).*A0(:,E.ind_C)), sum(P_IC(:,1).*A0(:,E.ind_H)), sum(P_IC(:,1).*A0(:,E.ind_O)), sum(P_IC(:,1).*A0(:,E.ind_N))])
        % Correction of the number of moles of CO2, H2O, O2 and N2 from atom
        % conservation equations
        NCO2_old = NCO2;
        NCO2 = NCO2_0-sum(P_IC(:,1).*A0(:,E.ind_C)); % C-atom conservation
        if TP > 3000, NCO2 = NCO2_old+relax*(NCO2-NCO2_old); end
        if NCO2 < 0, NCO2 = 0.75*NCO2_old; end
        
        NH2O_old = NH2O;
        NH2O = NH2O_0-sum(P_IC(:,1).*A0(:,E.ind_H))/2; % H-atom conservation
        if TP > 3000, NH2O = NH2O_old+relax*(NH2O-NH2O_old); end
        if NH2O < 0, NH2O = 0.75*NH2O_old; end
        
        P_IC(S.ind_CO2,1)= NCO2;
        P_IC(S.ind_H2O,1)= NH2O;
        
        NO2_old = NO2;
        NO2  = NO2_0+NCO2_0+NCO_0/2+NH2O_0/2-sum(P_IC(:,1).*A0(:,E.ind_O))/2; % O-atom conservation
        if TP > 3000, NO2 = NO2_old+relax*(NO2-NO2_old); end
        if NO2 < 0, NO2 = 0.75*NO2_old; end
        
        P_IC(S.ind_O2,1) = NO2;
        NN2_old = NN2;
        NN2  = NN2_0 -sum(P_IC(:,1).*A0(:,E.ind_N))/2; % N-atom conservation
        if TP > 3000, NN2 = NN2_old+relax*(NN2-NN2_old); end
        if NN2 < 0, NN2 = 0.75*NN2_old; end
        
        P_IC(S.ind_N2,1) = NN2;
        P_IC(S.ind_He,1) = NHe;
        P_IC(S.ind_Ar,1) = NAr;
        % Overall number of moles in the product species
        NP = sum(P_IC(:,1).*(1-P_IC(:,10)));
%         NP = sum(P_IC(:,1));
        DeltaNP = norm([NP-NP_old, x-sum(P_IC(:,1).*A0(:,E.ind_C)), y-sum(P_IC(:,1).*A0(:,E.ind_H)), z-sum(P_IC(:,1).*A0(:,E.ind_O)), w-sum(P_IC(:,1).*A0(:,E.ind_N))]);
        % relax = abs(relax*(1-DeltaNP/NP));
        % disp(sortrows([cell2table(S.NameSpecies(P_IC(:,1)>0.000001), 'VariableNames', {'Species'}), array2table(P_IC(P_IC(:,1)>0.000001), 'VariableNames', {'N'})], 2, 'descend'))
        end
    end % PHI <= 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [P_IC,DeltaNP] = incomplete_phi_2(P_IC,E,S,M,C)
        DeltaNP = 1;
        while abs(DeltaNP/NP) > C.tolN && it<itMax 
            it = it+1;
            % Initialization of the product matrix for incomplete combustion (IC)
            NP_old   = NP;
            P_IC_old = P_IC;
            
            P_IC = M0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % In rich combustion without carbon atoms the product mixture
            % contains H2O, H2 and N2 (in fuel-air combustion). The number
            % of moles of these species must be calculated using the H, O
            % and N-atom conservation equations
            %
            % For the minor species, we calculate the number of moles from
            % the equilibrium condition for the generalized reaction
            %
            % C.alpha CO2+(C.gamma-2*C.alpha) H2O+(C.beta/2-C.gamma+2*C.alpha) H2
            % +(C.omega/2) N2 <-VI-> C_alpha H_beta O_gamma N_omega [k6]
            %
            % with C.alpha = 0 if no carbon atoms are present.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Determination of the number of moles of the minor species from
            % the equilibrium condition for the above reaction
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if M.Lminors>0
                Ni     = k6.* NCO2.^C.alpha.* NH2O.^(C.gamma-2*C.alpha) ...
                    .* NH2.^(C.beta/2-C.gamma+2*C.alpha) .* NN2.^(C.omega/2) ...
                    .* zeta.^DNfactor_VI;
%                 Ni     = exp(log(k6)+log(NCO2).*(C.alpha)+log(NH2O).*(C.gamma-2*C.alpha) ...
%                     +log(NH2).*(C.beta/2-C.gamma(n)+2*C.alpha)+log(NN2+1).*(C.omega/2) ...
%                     +log(zeta).*DNfactor_VI);
%                 if M.major_OH && DeltaNP
% %                     Ni(M.ind_m_OH) = sqrt(NH2*NO2)/k11;
%                     Ni(M.ind_m_OH) = sqrt(NH2*NO2*k11*zeta^(-3/2));
%                 end
                Ni_old = P_IC_old(M.ind_minor,1)';
                aux = find(Ni_old~=0);
                if ~isempty(aux)
                    Ni(aux) = Ni_old(aux)+relax*(Ni(aux)-Ni_old(aux));
                end
                Ni(Ni>NP_old) = 0;
                P_IC(M.ind_minor,1) = Ni;
            end
            NO2     = zeta*(k2*NH2O/NH2)^2;
            NO2_old = P_IC_old(S.ind_O2,1);
            if NO2_old ~=0
                NO2 = NO2_old+relax*(NO2-NO2_old);
            end
            
            % Correction of the number of moles of H2, H2O and N2 from
            % atom conservation equations
            
            NH2O = NH2O_0-2*(NO2+sum(P_IC(:,1).*A0(:,E.ind_O))/2);
            NH2  = NH2_0 +2*(NO2+sum(P_IC(:,1).*A0(:,E.ind_O))/2-sum(P_IC(:,1).*A0(:,E.ind_H))/4);
            NN2  = NN2_0 -sum(P_IC(:,1).*A0(:,E.ind_N))/2; % N-atom conservation

            P_IC(S.ind_H2O,1)= NH2O;
            P_IC(S.ind_H2,1) = NH2;
            P_IC(S.ind_O2,1) = NO2;
            P_IC(S.ind_N2,1) = NN2;
            P_IC(S.ind_He,1) = NHe;
            P_IC(S.ind_Ar,1) = NAr;
            % Overall number of moles in the product species
            NP = sum(P_IC(:,1).*(1-P_IC(:,10)));
%             NP = sum(P_IC(:,1));
%             DeltaNP = NP-NP_old;
            DeltaNP = norm([NP-NP_old, y-sum(P_IC(:,1).*A0(:,E.ind_H)), z-sum(P_IC(:,1).*A0(:,E.ind_O)), w-sum(P_IC(:,1).*A0(:,E.ind_N))]);
            % relax = abs(relax*(1-DeltaNP/NP));
        end
    end % 1 < PHI < PHI_C AND ~H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [P_IC,DeltaNP] = incomplete_phi_3(P_IC,E,S,M,C)
        DeltaNP = 1;
        while abs(DeltaNP/NP) > C.tolN && it<itMax 
            it = it+1;
            % Initialization of the product matrix for incomplete combustion (IC)
            NP_old   = NP;
            P_IC_old = P_IC;
            
            P_IC = M0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % In rich combustion without hydrogen atoms the product mixture
            % contains CO2, CO and N2 (in fuel-air combustion). The number
            % of moles of these species must be calculated using the C, O
            % and N-atom conservation equations
            %
            % For the minor species, we calculate the number of moles from
            % the equilibrium condition for the generalized reaction
            %
            % (C.gamma-C.alpha-C.beta/2) CO2+(C.beta/2) H2O
            %    +(2*C.alpha-C.gamma+C.beta/2) CO+(C.omega/2) N2
            %                 <-V-> C_alpha H_beta O_gamma N_omega   [k5]
            %
            % with C.beta = 0 if no hydrogen atoms are present.
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Determination of the number of moles of the minor species from
            % the equilibrium condition for the above reaction
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if M.Lminors>0
                Ni     = k5 .* NCO2.^(C.gamma-C.alpha-C.beta/2) .* NH2O.^(C.beta/2) ...
                    .* NCO.^(C.beta/2-C.gamma+2*C.alpha) .* NN2.^(C.omega/2) ...
                    .* zeta.^DNfactor_V;
%                 Ni    = exp(log(k5)+(C.gamma-C.alpha-C.beta/2).*log(NCO2)+(C.beta/2).*log(NH2O) ...
%                        +log(NCO).*(2*C.alpha-C.gamma+C.beta/2)+log(NN2+1).*(C.omega/2) ...
%                        +log(zeta).*DNfactor_V);
                Ni_old = P_IC_old(M.ind_minor,1)';
                aux = find(Ni_old~=0);
                if ~isempty(aux)
                    Ni(aux) = Ni_old(aux)+relax*(Ni(aux)-Ni_old(aux));
                end
                Ni(Ni>NP_old) = 0;
                P_IC(M.ind_minor,1) = Ni;
            end
            NO2     = zeta*(k1*NCO2/NCO)^2;
            NO2_old = P_IC_old(S.ind_O2,1);
            if NO2_old ~=0
                NO2 = NO2_old+relax*(NO2-NO2_old);
            end
            
            % Correction of the number of moles of CO2, CO and N2 from
            % atom conservation equations
            
            NCO2 = NCO2_0-2*(NO2-3*sum(P_IC(:,1).*A0(:,E.ind_C))/2+sum(P_IC(:,1).*A0(:,E.ind_O))/2);
            NCO  = NCO_0 +2*(NO2-sum(P_IC(:,1).*A0(:,E.ind_C))+sum(P_IC(:,1).*A0(:,E.ind_O))/2);
            NN2  = NN2_0 -sum(P_IC(:,1).*A0(:,E.ind_N))/2; % N-atom conservation
            
            P_IC(S.ind_CO,1) = NCO;
            P_IC(S.ind_CO2,1)= NCO2;
            P_IC(S.ind_O2,1) = NO2;
            P_IC(S.ind_N2,1) = NN2;
            P_IC(S.ind_He,1) = NHe;
            P_IC(S.ind_Ar,1) = NAr;
            % Overall number of moles in the product species
            NP = sum(P_IC(:,1).*(1-P_IC(:,10)));
%             NP = sum(P_IC(:,1));
%             DeltaNP = NP-NP_old;
            DeltaNP = norm([NP-NP_old, x-sum(P_IC(:,1).*A0(:,E.ind_C)), z-sum(P_IC(:,1).*A0(:,E.ind_O)), w-sum(P_IC(:,1).*A0(:,E.ind_N))]);
            % relax = abs(relax*(1-DeltaNP/NP));
        end
    end % 1 < PHI < PHI_C AND ~C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [P_IC,DeltaNP] = incomplete_phi_4(P_IC,E,S,M,C)
                DeltaNP = 1;
        while abs(DeltaNP/NP) > C.tolN && it<itMax 
            it = it+1;
            % Initialization of the product matrix for incomplete combustion (IC)
            NP_old   = NP;
            P_IC_old = P_IC;
            
            P_IC = M0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % In rich combustion the product mixture always contains CO2, CO,
            % H2O, H2 and N2 (in fuel-air combustion). The number of moles of
            % these species must be calculated combining the 4 atom
            % conservation equations with the equilibrium condition for
            % the inverse water-gas shift reaction
            %
            %                  CO2+H2 <-IV-> H2O+CO             [k4]
            %
            % For the remaining minor species, we calculate the number of
            % moles from the equilibrium condition for the generalized
            % reactions
            %
            % (C.gamma-C.alpha-C.beta/2) CO2+(C.beta/2) H2O
            %    +(2*C.alpha-C.gamma+C.beta/2) CO+(C.omega/2) N2
            %                 <-V-> C_alpha H_beta O_gamma N_omega   [k5]
            %
            % or
            %
            % C.alpha CO2+(C.gamma-2*C.alpha) H2O+(C.beta/2-C.gamma+2*C.alpha) H2
            % +(C.omega/2) N2 <-VI-> C_alpha H_beta O_gamma N_omega [k6]
            %
            % Reaction V is preferred for low hydrogen content (e.g., CO
            % combustion), whereas reaction VI is preferred for low carbon
            % content (e.g., H2 combustion). In general, we shall use
            % reaction V, and leave reaction VI exclusively for H2
            % combustion.
            %
            % To estimate the ammount of O2 in the product mixture we use the
            % equilibrium condition for the reactions
            %
            % low Hydrogen:    CO2 <-I -> CO+(1/2) O2              [k1]
            % low Carbon:      H2O <-II-> H2+(1/2) O2             [k2]
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Determination of the number of moles of the minor species from
            % the equilibrium condition for the above reaction
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if M.Lminors>0 
                Ni     = k5 .* NCO2.^(C.gamma-C.alpha-C.beta/2) .* NH2O.^(C.beta/2) ...
                    .* NCO.^(2*C.alpha-C.gamma+C.beta/2) .* NN2.^(C.omega/2) ...
                    .* zeta.^DNfactor_V;
%                 Ni    = exp(log(k5)+(C.gamma-C.alpha-C.beta/2).*log(NCO2)+(C.beta/2).*log(NH2O) ...
%                        +log(NCO).*(2*C.alpha-C.gamma+C.beta/2)+log(NN2+1).*(C.omega/2) ...
%                        +log(zeta).*DNfactor_V);
                if M.major_CH4 && DeltaNP
%                         Ni(M.ind_m_CH4) = k8*NH2^2/zeta;
                        Ni(M.ind_m_CH4) = NH2*Ni(M.ind_m_CH3)/(Ni(M.ind_m_H)*k8);
%                         Ni(M.ind_m_CH4) = Ni(M.ind_m_CH4) + 0.995*(k8*NH2^2/zeta-Ni(M.ind_m_CH4));
%                         Ni(M.ind_m_C6H6)= k9*Ni(M.ind_m_C2H2)^3/zeta^2;
                end
%                   if M.minor_C
%                       Ni(M.ind_m_C)= sqrt(NCO^2/(NO2*zeta*k10));
%                   end
                Ni_old = P_IC_old(M.ind_minor,1)';
                aux = find(Ni_old~=0);
                if ~isempty(aux)
                    Ni(aux) = Ni_old(aux)+relax*(Ni(aux)-Ni_old(aux));
                end
                Ni(Ni>NP_old) = 0;
                P_IC(M.ind_minor,1) = Ni;
            end
            % Check
%             sum(P_IC(:,1).*A0(:,E.ind_C))
%             sum(P_IC(:,1).*A0(:,E.ind_H))
%             sum(P_IC(:,1).*A0(:,E.ind_O))
%             sum(P_IC(:,1).*A0(:,E.ind_N))
            %%%%
            NO2     = zeta*(k1*NCO2/NCO)^2;
            NO2_old = P_IC_old(S.ind_O2,1);
            if NO2_old ~=0
                NO2 = NO2_old+relax*(NO2-NO2_old);
            end
            
            % Correction of the number of moles of CO, H2, CO2, H2O and N2 from
            % atom conservation equations
            
            a = NCO2_0+NCO_0-sum(P_IC(:,1).*A0(:,E.ind_C));
            b = NH2O_0+NH2_0-sum(P_IC(:,1).*A0(:,E.ind_H))/2;
            c = NCO_0+NH2_0+2*(NO2-sum(P_IC(:,1).*A0(:,E.ind_C))-sum(P_IC(:,1).*A0(:,E.ind_H))/4+sum(P_IC(:,1).*A0(:,E.ind_O))/2);
            d = b-c;
            
            NCO2_old = NCO2;
            NCO_old  = NCO;
            NH2O_old = NH2O;
            NH2_old  = NH2;
            NN2_old  = NN2;
            
            NCO  = (1/2)*((a+c)*k4+d-sqrt((a-c)^2*k4^2+(2*(2*a*c+d*(a+c)))*k4+d^2))/(k4-1);
            NH2  = c-NCO;
            NH2O = d+NCO;
            NCO2 = a-NCO;
            NN2  = NN2_0 -sum(P_IC(:,1).*A0(:,E.ind_N))/2; % N-atom conservation
            
            if TP > 3000
                NCO2 = NCO2_old+0.1*relax*(NCO2-NCO2_old);
                NCO  = NCO_old+0.1*relax*(NCO -NCO_old);
                NH2O = NH2O_old+0.1*relax*(NH2O-NH2O_old);
                NH2  = NH2_old +0.1*relax*(NH2 -NH2_old);
                NN2  = NN2_old +0.1*relax*(NN2 -NN2_old);
            end
             
            if NCO2 < 0, NCO2 = 0.75*NCO2_old; end
            if NCO  < 0, NCO  = 0.75*NCO_old ; end
            if NH2O < 0, NH2O = 0.75*NH2O_old; end
            if NH2  < 0, NH2  = 0.75*NH2_old ; end
            if NN2  < 0, NN2  = 0.75*NN2_old ; end
            
            P_IC(S.ind_CO,1) = NCO;
            P_IC(S.ind_CO2,1)= NCO2;
            P_IC(S.ind_H2O,1)= NH2O;
            P_IC(S.ind_H2,1) = NH2;
            P_IC(S.ind_O2,1) = NO2;
            P_IC(S.ind_N2,1) = NN2;
            P_IC(S.ind_He,1) = NHe;
            P_IC(S.ind_Ar,1) = NAr;
            
            % Overall number of moles in the product species
            
            NP = sum(P_IC(:,1).*(1-P_IC(:,10)));
%             NP = sum(P_IC(:,1));
            DeltaNP = norm([NP-NP_old, x-sum(P_IC(:,1).*A0(:,E.ind_C)), y-sum(P_IC(:,1).*A0(:,E.ind_H)), z-sum(P_IC(:,1).*A0(:,E.ind_O)), w-sum(P_IC(:,1).*A0(:,E.ind_N))]);
%             if phi> 1.1, % relax = abs(relax*(1-DeltaNP/NP)); end
        end 
%         figure(1)
%         hold on;
%         T = linspace(300,6000);
%         r = 0.3385775+(0.9854897-0.3385775)./(1+(T./4058911).^1.817875).^658457.8;
%         plot(T,r);
%         pause(0.5);
    end % 1 < PHI < PHI_C OR ~SOOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [P_IC,DeltaNP] = incomplete_phi_5(P_IC,E,S,M,C)
        error('NOT IMPLEMENTED YET');
%         DeltaNP = 1;
%         while abs(DeltaNP/NP) > C.tolN
%             % Initialization of the product matrix for incomplete combustion (IC)
%             NP_old   = NP;
%             P_IC_old = P_IC;
%             
%             P_IC = M0;
%         end
    end % PHI <= 1 NOT IMPLEMENTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [P_IC,DeltaNP] = incomplete_phi_6(P_IC,E,S,M,C)
        DeltaNP = 1;
        while abs(DeltaNP/NP) > C.tolN && it<itMax 
            it = it+1;
            % Initialization of the product matrix for incomplete combustion (IC)
            NP_old   = NP;
            P_IC_old = P_IC;
                      
            P_IC = M0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % In rich combustion the product mixture always contains CO2, CO,
            % H2O, H2 and N2 (in fuel-air combustion). The number of moles of
            % these species must be calculated combining the 4 atom
            % conservation equations with the equilibrium condition for
            % the inverse water-gas shift reaction
            %
            %                  CO2+H2 <-IV-> H2O+CO             [k4]
            %
            % For the remaining minor species, we calculate the number of
            % moles from the equilibrium condition for the generalized
            % reactions
            %
            % (C.gamma-C.alpha-C.beta/2) CO2+(C.beta/2) H2O
            %    +(2*C.alpha-C.gamma+C.beta/2) CO+(C.omega/2) N2
            %                 <-V-> C_alpha H_beta O_gamma N_omega   [k5]
            %
            % or
            %
            % C.alpha CO2+(C.gamma-2*C.alpha) H2O+(C.beta/2-C.gamma+2*C.alpha) H2
            % +(C.omega/2) N2 <-VI-> C_alpha H_beta O_gamma N_omega [k6]
            %
            % Reaction V is preferred for low hydrogen content (e.g., CO
            % combustion), whereas reaction VI is preferred for low carbon
            % content (e.g., H2 combustion). In general, we shall use
            % reaction V, and leave reaction VI exclusively for H2
            % combustion.
            %
            % To estimate the ammount of O2 in the product mixture we use the
            % equilibrium condition for the reactions
            %
            % low Hydrogen:    CO2 <-I -> CO+(1/2) O2              [k1]
            % low Carbon:      H2O <-II-> H2+(1/2) O2             [k2]
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Determination of the number of moles of the minor species from
            % the equilibrium condition for the above reaction
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if t
%             if phi~=phi_c
%                 if M.Lminors>0
%                     Ni    = k8 .* NCO.^(C.gamma) .* NCgr.^(C.alpha-C.gamma) ...
%                         .* NH2.^(C.beta/2) .* NN2.^(C.omega/2) ...
%                         .* zeta.^DNfactor_VIII;
%                     for n = M.Lminors:-1:1
%                         Ni_old = P_IC_old(M.ind_minor(n),1);
%                         if Ni_old ~=0
%                             Ni(n) = Ni_old+relax*(Ni(n)-Ni_old);
%                         end
%                         P_IC(M.ind_minor(n),[strThProp.(minors_products{n}).Element_matrix(1,:),1]) = [Ni(n)*strThProp.(minors_products{n}).Element_matrix(2,:),Ni(n)];
%                     end
%                 end
%             else
                if M.Lminors>0
                    Ni    = k5 .* NCO2.^(C.gamma-C.alpha-C.beta/2) .* NH2O.^(C.beta/2) ...
                        .* NCO.^(2*C.alpha-C.gamma+C.beta/2) .* NN2.^(C.omega/2) ...
                        .* zeta.^DNfactor_V;
%                     Ni    = exp(log(k5)+(C.gamma-C.alpha-C.beta/2).*log(NCO2)+(C.beta/2).*log(NH2O) ...
%                        +log(NCO).*(2*C.alpha-C.gamma+C.beta/2)+log(NN2+1).*(C.omega/2) ...
%                        +log(zeta).*DNfactor_V);
                    if M.major_CH4 && DeltaNP
%                         Ni(M.ind_m_CH4) = k8*NH2^2/zeta;
                        Ni(M.ind_m_CH4) = NH2*Ni(M.ind_m_CH3)/(Ni(M.ind_m_H)*k8);
%                         Ni(M.ind_m_C) = Ni(M.ind_m_H)*Ni(M.ind_m_CH)/NH2*k10;
%                         Ni(M.ind_m_CH4) = Ni(M.ind_m_CH4) + 0.995*(k8*NH2^2/zeta-Ni(M.ind_m_CH4));
%                         Ni(M.ind_m_C6H6)= k9*Ni(M.ind_m_C2H2)^3/zeta^2;
                    end
                    Ni_old = P_IC_old(M.ind_minor,1)';
                    aux = find(Ni_old~=0);
                    if ~isempty(aux)
                        Ni(aux) = Ni_old(aux)+relax*(Ni(aux)-Ni_old(aux));
                    end
%                     Ni(Ni>NP_old) = 0;
                    P_IC(M.ind_minor,1) = Ni;
                end
%             end
%             NO2     = mean([zeta*(k1*NCO2/NCO)^2, zeta*(k2*NH2O/NH2)^2]);
            NO2     = zeta*(k1*NCO2/NCO)^2;
%             NO2 = 0;
            NO2_old = P_IC_old(S.ind_O2,1);
            if NO2_old ~= 0
                NO2 = NO2_old+relax*(NO2-NO2_old);
            end
            
            % Correction of the number of moles of CO, H2, CO2, H2O and N2 from
            % atom conservation equations
            
            a = NCO2_0+NCO_0+NCgr_0-sum(P_IC(:,1).*A0(:,E.ind_C));
            b = NH2O_0+NH2_0-sum(P_IC(:,1).*A0(:,E.ind_H))/2;
            c = 2*NCO2_0+NCO_0+NH2O_0-2*NO2-sum(P_IC(:,1).*A0(:,E.ind_O));
            
            %%% NEW
%             a0 = 4*a^3*k7^2+4*a^2*b*k7^2+a*b^2*k7^2-4*a^2*c*k7^2-2*a*b*c*k7^2+a*c^2*k7^2+(a*k4)/zeta^2-(c*k4)/zeta^2-(a*k4^2)/zeta^2+(c*k4^2)/zeta^2+(a^2*k7)/zeta+(a*b*k7)/zeta-(a*c*k7)/zeta+(4*a^2*k4*k7)/zeta-(4*a*c*k4*k7)/zeta-(b*c*k4*k7)/zeta+(c^2*k4*k7)/zeta-(4*a^2*k4^2*k7)/zeta+(4*a*c*k4^2*k7)/zeta-(c^2*k4^2*k7)/zeta;
%             a1 = 12*a^2*k7^2+8*a*b*k7^2+b^2*k7^2-8*a*c*k7^2-2*b*c*k7^2+c^2*k7^2+k4/zeta^2-k4^2/zeta^2+(2*a*k7)/zeta+(b*k7)/zeta-(c*k7)/zeta+(8*a*k4*k7)/zeta-(4*c*k4*k7)/zeta-(8*a*k4^2*k7)/zeta+(4*c*k4^2*k7)/zeta;
%             a2 = (-12*a*k7-4*b*k7+4*c*k7-1/zeta-(4*k4)/zeta+(4*k4^2)/zeta)*k7;
%             a3 = 4*k7^2;
            %%% NEW
            
            NCO2_old = NCO2;
            NCO_old  = NCO;
            NCgr_old = NCgr;
            NH2O_old = NH2O;
            NH2_old  = NH2;
            NN2_old  = NN2;
                       
            NCgr = real(1/24*(24*a+(4*(2+k4))/(k4*mu)+(2*(1+1i*sqrt(3))*(-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu)))/(k4*(8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))+(2*1i*(1i+sqrt(3))*(8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))/(k4*mu^2)-(1/(6*k4^2*mu^3))*((2*(2+k4)*mu+((1+1i*sqrt(3))*mu^2*(-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu)))/(8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3)+1i*(1i+sqrt(3))*(8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))^2)));  
%             NCgr = real(1/24*(24*a+(4*(2+k4))/(k4*mu)+(2*(1-1i*sqrt(3))*(-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu)))/(k4*(8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))-(2*1i*(-1i+sqrt(3))*(8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))/(k4*mu^2)-(1/(6*k4^2*mu^3))*((-2*(2+k4)*mu+(1i*(1i+sqrt(3))*mu^2*(-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu)))/(8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3)+(1+1i*sqrt(3))*(8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))^2)));  
            %%% NEW
%             NCgr = real(-(a2/(3*a3))+(2^(1/3)*(-a2^2+3*a1*a3))/(3*a3*(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3))-(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3)/(3*2^(1/3)*a3));
%             NCgr = real(-(a2/(3*a3))-((1-1i*sqrt(3))*(-a2^2+3*a1*a3))/(3*2^(2/3)*a3*(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3))+((1+1i*sqrt(3))*(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3))/(6*2^(1/3)*a3));
%             NCgr = real(-(a2/(3*a3))-((1+1i*sqrt(3))*(-a2^2+3*a1*a3))/(3*2^(2/3)*a3*(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3))+((1-1i*sqrt(3))*(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3))/(6*2^(1/3)*a3));
            %%% NEW
            NCO2 = (1+2*a*mu-2*mu*NCgr-sqrt(1+4*a*mu-4*mu*NCgr))/(2*mu);
            %%% NEW 2 
%             a0 = -2*k4/zeta-k7*b+2*k7*c;
%             a1 = 2*k7+4*k4*k7;
%             a2 = 4*k7^2*zeta;
%             a3 = 2*k4*c/zeta;
%             NCO = real(-(a1/(3*a2))-(2^(1/3)*(-a1^2-3*a0*a2))/(3*a2*(-2*a1^3-9*a0*a1*a2+27*a2^2*a3+sqrt(-4*(a1^2+3*a0*a2)^3+(2*a1^3+9*a0*a1*a2-27*a2^2*a3)^2))^(1/3))+(-2*a1^3-9*a0*a1*a2+27*a2^2*a3+sqrt(-4*(a1^2+3*a0*a2)^3+(2*a1^3+9*a0*a1*a2-27*a2^2*a3)^2))^(1/3)/(3*2^(1/3)*a2));
%             NCO2 = NCO^2*mu;
            %%% NEW 2
            if ~isreal(NCO2)  
            	NCgr = real(1/6*(6*a+(2+k4)/(k4*mu)+(4-2*k4+k4^2*(1-6*b*mu+6*c*mu))/(k4*(8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))+(1/(k4*mu^2))*((8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))-(1/(6*k4^2*mu^3))*((2*mu+k4*mu-(mu^2*(-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu)))/(8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3)+(8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))^2)));
                %%% NEW
%                 NCgr = real(-(a2/(3*a3))+(2^(1/3)*(-a2^2+3*a1*a3))/(3*a3*(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3))-(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3)/(3*2^(1/3)*a3));
                %%% NEW
                NCO2 = (1+2*a*mu-2*mu*NCgr-sqrt(1+4*a*mu-4*mu*NCgr))/(2*mu);
            end
            NCO = (-1+sqrt(1+4*a*mu-4*mu*NCgr))/(2*mu);
            NH2O = c-2*NCO2-NCO;
            NH2 = b-NH2O;
            %%%%%%%%%%%%%%%%%%%%%%%%
%             NH2O = real((1/(6*k4))*(2*(k4*(2*b+c-2*mu2)+mu2)+(2*2^(1/3)*(b^2*k4^2+(c*k4+mu2-2*k4*mu2)^2+b*k4*(mu2-2*k4*(c+4*mu2))))/(-2*b^3*k4^3+3*b^2*k4^2*(2*k4*(c-5*mu2)-mu2)+2*(c*k4+mu2-2*k4*mu2)^3-3*b*k4*(c*k4+mu2-2*k4*mu2)*(2*c*k4+(-1+8*k4)*mu2)+3*sqrt(3)*sqrt(-b^2*k4^2*mu2*(-8*b^3*k4^4+(8*c*k4^2+mu2)*(c*k4+mu2-2*k4*mu2)^2+2*b*k4*(-12*c^2*k4^3+c*k4*(-1+2*k4-40*k4^2)*mu2+(1+4*(-2+k4)*k4)*mu2^2)+b^2*k4^2*(mu2+4*k4*(6*c*k4+(-5+k4)*mu2)))))^(1/3)+2^(2/3)*(-2*b^3*k4^3+3*b^2*k4^2*(2*k4*(c-5*mu2)-mu2)+2*(c*k4+mu2-2*k4*mu2)^3-3*b*k4*(c*k4+mu2-2*k4*mu2)*(2*c*k4+(-1+8*k4)*mu2)+3*sqrt(3)*sqrt(-b^2*k4^2*mu2*(-8*b^3*k4^4+(8*c*k4^2+mu2)*(c*k4+mu2-2*k4*mu2)^2+2*b*k4*(-12*c^2*k4^3+c*k4*(-1+2*k4-40*k4^2)*mu2+(1+4*(-2+k4)*k4)*mu2^2)+b^2*k4^2*(mu2+4*k4*(6*c*k4+(-5+k4)*mu2)))))^(1/3)));
% %             NH2O = real((1/(12*k4))*(4*(k4*(2*b+c-2*mu2)+mu2)-(2*(-1i)*2^(1/3)*(-(-1i)+sqrt(3))*(b^2*k4^2+(c*k4+mu2-2*k4*mu2)^2+b*k4*(mu2-2*k4*(c+4*mu2))))/(-2*b^3*k4^3+3*b^2*k4^2*(2*k4*(c-5*mu2)-mu2)+2*(c*k4+mu2-2*k4*mu2)^3-3*b*k4*(c*k4+mu2-2*k4*mu2)*(2*c*k4+(-1+8*k4)*mu2)+3*sqrt(3)*sqrt(-b^2*k4^2*mu2*(-8*b^3*k4^4+(8*c*k4^2+mu2)*(c*k4+mu2-2*k4*mu2)^2+2*b*k4*(-12*c^2*k4^3+c*k4*(-1+2*k4-40*k4^2)*mu2+(1+4*(-2+k4)*k4)*mu2^2)+b^2*k4^2*(mu2+4*k4*(6*c*k4+(-5+k4)*mu2)))))^(1/3)+(-1i)*2^(2/3)*((-1i)+sqrt(3))*(-2*b^3*k4^3+3*b^2*k4^2*(2*k4*(c-5*mu2)-mu2)+2*(c*k4+mu2-2*k4*mu2)^3-3*b*k4*(c*k4+mu2-2*k4*mu2)*(2*c*k4+(-1+8*k4)*mu2)+3*sqrt(3)*sqrt(-b^2*k4^2*mu2*(-8*b^3*k4^4+(8*c*k4^2+mu2)*(c*k4+mu2-2*k4*mu2)^2+2*b*k4*(-12*c^2*k4^3+c*k4*(-1+2*k4-40*k4^2)*mu2+(1+4*(-2+k4)*k4)*mu2^2)+b^2*k4^2*(mu2+4*k4*(6*c*k4+(-5+k4)*mu2)))))^(1/3)));
%             NH2 = b-NH2O;
%             if NH2 < 0 || ~isreal(NH2)
%                 NH2O = real((1/(12*k4))*(4*(k4*(2*b+c-2*mu2)+mu2)-(2*1i*2^(1/3)*(-1i+sqrt(3))*(b^2*k4^2+(c*k4+mu2-2*k4*mu2)^2+b*k4*(mu2-2*k4*(c+4*mu2))))/(-2*b^3*k4^3+3*b^2*k4^2*(2*k4*(c-5*mu2)-mu2)+2*(c*k4+mu2-2*k4*mu2)^3-3*b*k4*(c*k4+mu2-2*k4*mu2)*(2*c*k4+(-1+8*k4)*mu2)+3*sqrt(3)*sqrt(-b^2*k4^2*mu2*(-8*b^3*k4^4+(8*c*k4^2+mu2)*(c*k4+mu2-2*k4*mu2)^2+2*b*k4*(-12*c^2*k4^3+c*k4*(-1+2*k4-40*k4^2)*mu2+(1+4*(-2+k4)*k4)*mu2^2)+b^2*k4^2*(mu2+4*k4*(6*c*k4+(-5+k4)*mu2)))))^(1/3)+1i*2^(2/3)*(1i+sqrt(3))*(-2*b^3*k4^3+3*b^2*k4^2*(2*k4*(c-5*mu2)-mu2)+2*(c*k4+mu2-2*k4*mu2)^3-3*b*k4*(c*k4+mu2-2*k4*mu2)*(2*c*k4+(-1+8*k4)*mu2)+3*sqrt(3)*sqrt(-b^2*k4^2*mu2*(-8*b^3*k4^4+(8*c*k4^2+mu2)*(c*k4+mu2-2*k4*mu2)^2+2*b*k4*(-12*c^2*k4^3+c*k4*(-1+2*k4-40*k4^2)*mu2+(1+4*(-2+k4)*k4)*mu2^2)+b^2*k4^2*(mu2+4*k4*(6*c*k4+(-5+k4)*mu2)))))^(1/3)));
%                 NH2 = b-NH2O;
%             end
%             NCO2 = mu2*(NH2O/NH2)^2;
%             NCO = c-2*NCO2-NH2O;
%             NCgr = a-NCO2-NCO;
            %%%%%%%%%%%%%%%%%%%%%%%%
            NN2  = NN2_0-sum(P_IC(:,1).*A0(:,E.ind_N))/2; % N-atom conservation
            
            if TP > 3000
                NCO2 = NCO2_old+relax*(NCO2-NCO2_old);
                NCO  = NCO_old +relax*(NCO -NCO_old);
                NCgr = NCgr_old+relax*(NCgr-NCgr_old);
                NH2O = NH2O_old+relax*(NH2O-NH2O_old);
                NH2  = NH2_old +relax*(NH2 -NH2_old);
                NN2  = NN2_old +relax*(NN2 -NN2_old);
            end
            
            end
%             if NCgr <= 0 || ~t
            if NCgr <= 1e-5   
                if t
                    NCgr = 0;

%                     NCO_0 = (-1+sqrt(1+4*a*mu-4*mu*NCgr))/(2*mu);
%                     NCO2_0 = (1+2*a*mu-2*mu*NCgr-sqrt(1+4*a*mu-4*mu*NCgr))/(2*mu);
%                     NH2O_0 = c-2*NCO2-NCO;
%                     NH2_0 = b-NH2O;

                    if it>1
                        NCO_0  = NCO_old;
                        NCO2_0 = NCO2_old;
                        NH2O_0 = NH2O_old;
                        NH2_0  = NH2_old;
                        NN2_0  = NN2_old;
                        NO2_0  = NO2_old;
                        NCgr_0 = NCgr_old;
                    
                        N_minor_C_0 = sum(P_IC_old(:,1).*A0(:,E.ind_C))-NCO2_0-NCO_0-NCgr_0;
                        N_minor_H_0 = sum(P_IC_old(:,1).*A0(:,E.ind_H))-2*NH2O_0-2*NH2_0;
                        N_minor_O_0 = sum(P_IC_old(:,1).*A0(:,E.ind_O))-2*NO2_0-2*NCO2_0-NCO_0-NH2O_0;
                        N_minor_N_0 = sum(P_IC_old(:,1).*A0(:,E.ind_N))-2*NN2_0;
                        
                        N_minor_C_0_nswt = sum(P_IC_old(:,1).*A0(:,E.ind_C).*(1-P_IC(:,10)))-NCO2_0-NCO_0-NCgr_0;
                        N_minor_H_0_nswt = sum(P_IC_old(:,1).*A0(:,E.ind_H).*(1-P_IC(:,10)))-2*NH2O_0-2*NH2_0;
                        N_minor_O_0_nswt = sum(P_IC_old(:,1).*A0(:,E.ind_O).*(1-P_IC(:,10)))-2*NO2_0-2*NCO2_0-NCO_0-NH2O_0;
                        N_minor_N_0_nswt = sum(P_IC_old(:,1).*A0(:,E.ind_N).*(1-P_IC(:,10)))-2*NN2_0;
                    else
%                         NCO_0  = NCO;
%                         NCO2_0 = NCO2;
%                         NH2O_0 = NH2O;
%                         NH2_0  = NH2;
%                         NN2_0  = NN2;
%                         NO2_0  = NO2;
%                         NCgr_0 = NCgr;
                        
%                         N_minor_C_0 = sum(P_IC(:,1).*A0(:,E.ind_C));
%                         N_minor_H_0 = sum(P_IC(:,1).*A0(:,E.ind_H));
%                         N_minor_O_0 = sum(P_IC(:,1).*A0(:,E.ind_O));
%                         N_minor_N_0 = sum(P_IC(:,1).*A0(:,E.ind_N));
%                         
%                         N_minor_C_0_nswt = sum(P_IC(:,1).*A0(:,E.ind_C).*(1-P_IC(:,10)));
%                         N_minor_H_0_nswt = sum(P_IC(:,1).*A0(:,E.ind_H).*(1-P_IC(:,10)));
%                         N_minor_O_0_nswt = sum(P_IC(:,1).*A0(:,E.ind_O).*(1-P_IC(:,10)));
%                         N_minor_N_0_nswt = sum(P_IC(:,1).*A0(:,E.ind_N).*(1-P_IC(:,10)));
                        
                        N_minor_C_0 = 0;
                        N_minor_H_0 = 0;
                        N_minor_O_0 = 0;
                        N_minor_N_0 = 0;

                        N_minor_C_0_nswt = 0;
                        N_minor_H_0_nswt = 0;
                        N_minor_O_0_nswt = 0;
                        N_minor_N_0_nswt = 0;
%                         NP_old = N_minor_C_0+N_minor_H_0+N_minor_O_0+N_minor_N_0...
%                             +NCO2_0+NCO_0+NH2O_0+NH2_0+NO2_0+NN2_0+NCgr_0...
%                             +NHe+NAr;
%                         NP_old = N_minor_C_0_nswt+N_minor_H_0_nswt+N_minor_O_0_nswt+N_minor_N_0_nswt...
%                             +NCO2_0+NCO_0+NH2O_0+NH2_0+NO2_0+NN2_0...
%                             +NHe+NAr;
                    end
                    
                    NP_old = N_minor_C_0_nswt+N_minor_H_0_nswt+N_minor_O_0_nswt+N_minor_N_0_nswt...
                            +y/4+(NH2_0+NCO_0+z+w)/2+NHe+NAr;
                    NP_old = 0.5*(NP+NP_old); 
                    
                    if strfind(PD.ProblemType,'P') == 2 % PD.ProblemType = 'TP', 'HP', or 'SP'
                        zeta = NP_old/pP;
                    elseif strfind(PD.ProblemType,'V') == 2 % PD.ProblemType = 'TV', 'EV', or 'SV'
                        zeta = (vP/1000)*1e5/R0TP; % vP is given in l
                    end
                    
                    NCO2 = NCO2_0;
                    NCO = NCO_0;
                    NH2O = NH2O_0;
                    NH2 = NH2_0;
                end
                
                %%%%%%%%%%%%%%%%
                % CHECK
%                 N_minor_C_0+NCO2_0+NCO_0+NCgr_0
%                 N_minor_H_0+2*NH2O_0+2*NH2_0 
%                 N_minor_O_0+2*NCO2_0+NCO_0+2*NO2_0+NH2O_0
%                 N_minor_N_0+2*NN2_0 
                %%%%%%%%%%%%%%%%
                if M.Lminors>0
                    Ni    = k5 .* NCO2.^(C.gamma-C.alpha-C.beta/2) .* NH2O.^(C.beta/2) ...
                        .* NCO.^(2*C.alpha-C.gamma+C.beta/2) .* NN2.^(C.omega/2) ...
                        .* zeta.^DNfactor_V;
%                     Ni    = exp(log(k5)+(C.gamma-C.alpha-C.beta/2).*log(NCO2)+(C.beta/2).*log(NH2O) ...
%                        +log(NCO).*(2*C.alpha-C.gamma+C.beta/2)+log(NN2).*(C.omega/2) ...
%                        +log(zeta).*DNfactor_V);
                    if M.major_CH4 && DeltaNP
%                         Ni(M.ind_m_CH4) = k8*NH2^2/zeta;
                    Ni(M.ind_m_CH4) = NH2*Ni(M.ind_m_CH3)/(Ni(M.ind_m_H)*k8);
%                         Ni(M.ind_m_CH4) = Ni(M.ind_m_CH4) + 0.995*(k8*NH2^2/zeta-Ni(M.ind_m_CH4));
%                         Ni(M.ind_m_C6H6)= k9*Ni(M.ind_m_C2H2)^3/zeta^2;
                    end
                    Ni_old = P_IC_old(M.ind_minor,1)';
                    aux = find(Ni_old~=0);
                    if ~isempty(aux)
                        Ni(aux) = Ni_old(aux)+relax*(Ni(aux)-Ni_old(aux));
                        %                           Ni(aux) = exp(log(Ni_old(aux)) +relax*(log(Ni(aux)) -log(Ni_old(aux))));
                    end
%                     Ni(Ni>NP_old) = 0;
                    P_IC(M.ind_minor,1) = Ni;
                end
                
                NO2     = zeta*(k1*NCO2/NCO)^2;
                
                NO2_old = P_IC_old(S.ind_O2,1);
                if NO2_old ~= 0
                    NO2 = NO2_old+relax*(NO2-NO2_old);
%                     NO2 = exp(log(NO2_old) +relax*(log(NO2) -log(NO2_old)));
                end
                
                a = NCO2_0+NCO_0+NCgr_0+N_minor_C_0-sum(P_IC(:,1).*A0(:,E.ind_C));
                b = NH2O_0+NH2_0+N_minor_H_0/2-sum(P_IC(:,1).*A0(:,E.ind_H))/2;
                c = NCO_0+NH2_0+2*(NCgr_0+N_minor_C_0+N_minor_H_0/4-N_minor_O_0/2)+2*(NO2-NO2_0-sum(P_IC(:,1).*A0(:,E.ind_C))-sum(P_IC(:,1).*A0(:,E.ind_H))/4+sum(P_IC(:,1).*A0(:,E.ind_O))/2);

                d = b-c;
                
                NCO2_old = NCO2;
                NCO_old  = NCO;
                NH2O_old = NH2O;
                NH2_old  = NH2;
                NN2_old  = NN2;
                
                NCO  = (1/2)*((a+c)*k4+d-sqrt((a-c)^2*k4^2+(2*(2*a*c+d*(a+c)))*k4+d^2))/(k4-1);
                NCO2 = a-NCO;
                NH2  = c-NCO;
                NH2O = d+NCO;

                NN2  = NN2_0+N_minor_N_0/2-sum(P_IC(:,1).*A0(:,E.ind_N))/2; % N-atom conservation
                
%                 NCO2 = NCO2_old+0.0001*relax*(NCO2-NCO2_old);
%                 NCO  = NCO_old +relax*(NCO -NCO_old);
%                 NH2O = NH2O_old+0.0001*relax*(NH2O-NH2O_old);
%                 NH2  = NH2_old +relax*(NH2 -NH2_old);

                NCO2 = NCO2_old+0.001*relax*(NCO2-NCO2_old);
%                 NCO  = NCO_old +0.01*relax*(NCO -NCO_old);
                NH2O = NH2O_old+0.001*relax*(NH2O-NH2O_old);
%                 NH2  = NH2_old +0.01*relax*(NH2 -NH2_old);
                
%                 NCO2 = NCO2_old+0.01*relax*(NCO2-NCO2_old);
%                 NCO  = NCO_old +0.01*relax*(NCO -NCO_old);
%                 NH2O = NH2O_old+0.01*relax*(NH2O-NH2O_old);
%                 NH2  = NH2_old +0.01*relax*(NH2 -NH2_old);

%                 NCO  = (NCO_old-NCO)/log(NCO_old/NCO);
%                 NCO2  = (NCO2_old-NCO2)/log(NCO2_old/NCO2);
%                 NH2O  = (NH2O_old-NH2O)/log(NH2O_old/NH2O);
%                 NH2  = (NH2_old-NH2)/log(NH2_old/NH2);
                
                if TP > 3000
                    NCO2 = NCO2_old+relax*(NCO2-NCO2_old);
                    NCO  = NCO_old +relax*(NCO -NCO_old);
                    NH2O = NH2O_old+relax*(NH2O-NH2O_old);
                    NH2  = NH2_old +relax*(NH2 -NH2_old);
                    NN2  = NN2_old +relax*(NN2 -NN2_old);
                    NCgr  = NCgr_old +relax*(NCgr -NCgr_old);
                end

                t = false;
            end
           
            
            if NCO2 < 0 || ~isreal(NCO2), NCO2 = 0.75*NCO2_old; end
            if NCO  < 0 || ~isreal(NCO), NCO  = 0.75*NCO_old ; end
            if NCgr < 0 || ~isreal(NCgr), NCgr = 0.01*NCgr_old ; end
            if NH2O < 0 || ~isreal(NH2O), NH2O = 0.75*NH2O_old; end
            if NH2  < 0 || ~isreal(NH2), NH2  = 0.75*NH2_old ; end
            if NN2  < 0 || ~isreal(NN2), NN2  = 0.75*NN2_old ; end
            
            P_IC(S.ind_Cgr,1)= NCgr;
            P_IC(S.ind_CO,1) = NCO;
            P_IC(S.ind_CO2,1)= NCO2;
            P_IC(S.ind_H2O,1)= NH2O;
            P_IC(S.ind_H2,1) = NH2;
            P_IC(S.ind_O2,1) = NO2;
            P_IC(S.ind_N2,1) = NN2;
            P_IC(S.ind_He,1) = NHe;
            P_IC(S.ind_Ar,1) = NAr;

            % Overall number of moles in the product species
            NP = sum(P_IC(:,1).*(1-P_IC(:,10)));
%             NP = sum(P_IC(:,1));
            DeltaNP = norm([NP-NP_old, x-sum(P_IC(:,1).*A0(:,E.ind_C)), y-sum(P_IC(:,1).*A0(:,E.ind_H)), z-sum(P_IC(:,1).*A0(:,E.ind_O)), w-sum(P_IC(:,1).*A0(:,E.ind_N))]);
            % relax = abs(relax*(1-DeltaNP/NP));
        end
    end % PHI >= PHI_C OR SOOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [P_IC,DeltaNP] = incomplete_phi_7(P_IC,E,S,M,C)
        DeltaNP = 1;
        while abs(DeltaNP/NP) > C.tolN && it<itMax 
            it = it+1;
            % Initialization of the product matrix for incomplete combustion (IC)
            NP_old   = NP;
            P_IC_old = P_IC;
                      
            P_IC = M0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % In rich combustion the product mixture always contains CO2, CO,
            % H2O, H2 and N2 (in fuel-air combustion). The number of moles of
            % these species must be calculated combining the 4 atom
            % conservation equations with the equilibrium condition for
            % the inverse water-gas shift reaction
            %
            %                  CO2+H2 <-IV-> H2O+CO             [k4]
            %
            % For the remaining minor species, we calculate the number of
            % moles from the equilibrium condition for the generalized
            % reactions
            %
            % (C.gamma-C.alpha-C.beta/2) CO2+(C.beta/2) H2O
            %    +(2*C.alpha-C.gamma+C.beta/2) CO+(C.omega/2) N2
            %                 <-V-> C_alpha H_beta O_gamma N_omega   [k5]
            %
            % or
            %
            % C.alpha CO2+(C.gamma-2*C.alpha) H2O+(C.beta/2-C.gamma+2*C.alpha) H2
            % +(C.omega/2) N2 <-VI-> C_alpha H_beta O_gamma N_omega [k6]
            %
            % Reaction V is preferred for low hydrogen content (e.g., CO
            % combustion), whereas reaction VI is preferred for low carbon
            % content (e.g., H2 combustion). In general, we shall use
            % reaction V, and leave reaction VI exclusively for H2
            % combustion.
            %
            % To estimate the ammount of O2 in the product mixture we use the
            % equilibrium condition for the reactions
            %
            % low Hydrogen:    CO2 <-I -> CO+(1/2) O2              [k1]
            % low Carbon:      H2O <-II-> H2+(1/2) O2             [k2]
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Determination of the number of moles of the minor species from
            % the equilibrium condition for the above reaction
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if t
                if M.Lminors>0
                    Ni    = k5 .* NCO2.^(C.gamma-C.alpha-C.beta/2) .* NH2O.^(C.beta/2) ...
                        .* NCO.^(2*C.alpha-C.gamma+C.beta/2) .* NN2.^(C.omega/2) ...
                        .* zeta.^DNfactor_V;
                    Ni(Ni>NP_old) = 0;
                    Ni_old = P_IC_old(M.ind_minor,1)';
                    aux = find(Ni_old~=0);
                    if ~isempty(aux)
                        Ni(aux) = Ni_old(aux)+relax*(Ni(aux)-Ni_old(aux));
                    end
                    P_IC(M.ind_minor,1) = Ni;
                end
            
            NO2     = zeta*(k1*NCO2/NCO)^2;
            NO2_old = P_IC_old(S.ind_O2,1);
            if NO2_old ~= 0
                NO2 = NO2_old+relax*(NO2-NO2_old);
            end
            
            % Correction of the number of moles of CO, H2, CO2, H2O and N2 from
            % atom conservation equations
            
            a = NCO2_0+NCO_0+NCgr_0-sum(P_IC(:,1).*A0(:,E.ind_C));
            b = NH2O_0+NH2_0-sum(P_IC(:,1).*A0(:,E.ind_H))/2;
            c = 2*NCO2_0+NCO_0+NH2O_0+2*NO2_0-2*NO2-sum(P_IC(:,1).*A0(:,E.ind_O));
            
            NCO2_old = NCO2;
            NCO_old  = NCO;
            NCgr_old = NCgr;
            NH2O_old = NH2O;
            NH2_old  = NH2;
            NN2_old  = NN2;
            
            mu = k7/zeta;
            
            NCgr = real(1/24*(24*a+(4*(2+k4))/(k4*mu)+(2*(1+1i*sqrt(3))*(-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu)))/(k4*(8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))+(2*1i*(1i+sqrt(3))*(8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))/(k4*mu^2)-(1/(6*k4^2*mu^3))*((2*(2+k4)*mu+((1+1i*sqrt(3))*mu^2*(-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu)))/(8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3)+1i*(1i+sqrt(3))*(8*mu^3-6*k4*mu^3-3*k4^2*mu^3+k4^3*mu^3-18*b*k4^2*mu^4-36*c*k4^2*mu^4-9*b*k4^3*mu^4+9*c*k4^3*mu^4+sqrt(mu^6*((-4+2*k4+k4^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*k4+k4^3*(-1+9*b*mu-9*c*mu)+3*k4^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))^2)));  
            NCO2 = (1+2*a*mu-2*mu*NCgr-sqrt(1+4*a*mu-4*mu*NCgr))/(2*mu);
            NCO = (-1+sqrt(1+4*a*mu-4*mu*NCgr))/(2*mu);
            NH2O = c-2*NCO2-NCO;
            NH2 = b-NH2O;
            
            NN2  = NN2_0-sum(P_IC(:,1).*A0(:,E.ind_N))/2; % N-atom conservation
            end
            
            
            if TP > 3000
                NCO2 = NCO2_old+relax*(NCO2-NCO2_old);
                NCO  = NCO_old +relax*(NCO -NCO_old);
                NCgr = NCgr_old+relax*(NCgr-NCgr_old);
                NH2O = NH2O_old+relax*(NH2O-NH2O_old);
                NH2  = NH2_old +relax*(NH2 -NH2_old);
                NN2  = NN2_old +relax*(NN2 -NN2_old);
            end
            
            if NCO2 < 0, NCO2 = 0.75*NCO2_old; end
            if NCO  < 0, NCO  = 0.75*NCO_old ; end
            if NCgr < 0, NCgr = 0.01*NCgr_old ; end
            if NH2O < 0, NH2O = 0.75*NH2O_old; end
            if NH2  < 0, NH2  = 0.75*NH2_old ; end
            if NN2  < 0, NN2  = 0.75*NN2_old ; end
            
            P_IC(S.ind_Cgr,1)= NCgr;
            P_IC(S.ind_CO,1) = NCO;
            P_IC(S.ind_CO2,1)= NCO2;
            P_IC(S.ind_H2O,1)= NH2O;
            P_IC(S.ind_H2,1) = NH2;
            P_IC(S.ind_O2,1) = NO2;
            P_IC(S.ind_N2,1) = NN2;
            P_IC(S.ind_He,1) = NHe;
            P_IC(S.ind_Ar,1) = NAr;
            
            % Overall number of moles in the product species
            
            NP = sum(P_IC(:,1).*(1-P_IC(:,10)));
%             NP = sum(P_IC(:,1));
            
            DeltaNP = norm([NP-NP_old, x-sum(P_IC(:,1).*A0(:,E.ind_C)), y-sum(P_IC(:,1).*A0(:,E.ind_H)), z-sum(P_IC(:,1).*A0(:,E.ind_O)), w-sum(P_IC(:,1).*A0(:,E.ind_N))]);
            
            if NCgr <= 1e-8
               P_IC = incomplete_phi_6(P_IC); 
               break; 
            end
        end
    end % TEST
end
