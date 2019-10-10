function N_IC = CalculateProductsIC(P_CC,TP,pP,vP,phi,minor_products,phi_c,FLAG_SOOT,M0,A0,E,S,C,M,PD,TN,strThProp)

R0TP = C.R0*TP;
it = 0; itMax = 500; t=true;
% Relaxation/iteration parameters
relax = 0.00007385775+(0.9854897-0.00007385775)/(1+(TP/4058911)^1.817875)^658457.8;
% Number of moles of the major species in the product mixture under the
% assumption of complete combustion (CC), denoted by subscript _0
NCO2_0 = P_CC(S.idx_CO2,1);
NCO_0  = P_CC(S.idx_CO,1);
NH2O_0 = P_CC(S.idx_H2O,1);
NH2_0  = P_CC(S.idx_H2,1);
NO2_0  = P_CC(S.idx_O2,1);
NN2_0  = P_CC(S.idx_N2,1);
NCgr_0 = P_CC(S.idx_Cgr,1);
% Number of C, H, O, N, He, Ar-atoms in the product species
NatomE = sum(P_CC(:,1).*A0);

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

NP = sum(P_CC(:,1).*(1-P_CC(:,10))); % Sum of num of moles of gases-(1-swt), with swt == condensed phase.
% NP = sum(P_CC(:,1));

if strfind(PD.ProblemType,'P') == 2 % PD.ProblemType = 'TP', 'HP', or 'SP'
    DNfactor = NP/pP;
elseif strfind(PD.ProblemType,'V') == 2 % PD.ProblemType = 'TV', 'EV', or 'SV'
    DNfactor = (vP/1000)*1e5/R0TP; % vP is given in l
end

N_IC = M0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATION OF GIBBS FREE ENERGY, CONSTANTS EQUILIBRIUM AND OTHERS
g_CO2 = species_g0_new('CO2',TP,strThProp);
g_CO  = species_g0_new('CO',TP,strThProp);
g_H2O = species_g0_new('H2O',TP,strThProp);
DG0_I    = (g_CO-g_CO2)*1000;
DG0_II   = -g_H2O*1000;
KPT_I  = exp(-DG0_I/R0TP);
KPT_II  = exp(-DG0_II/R0TP);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correction to the initial guess of oxygene moles and overall number of
% moles
NO2    = NO2_0+((NH2O*KPT_II+NCO2*KPT_I)/2)^(2/3)*DNfactor^(1/3);
NP     = NP+(NO2-NO2_0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if phi<=1 && M.L_minor>0 % case of lean-to-stoichiometric mixtures
    DNfactor_III = 1-(C.beta+2*(C.gamma+C.omega))/4;
    for n = M.L_minor:-1:1
        DG0_III(n) = (species_g0_new(minor_products{n},TP,strThProp)-C.alpha(n)*g_CO2 ...
            -(C.beta(n)/2)*g_H2O)*1000;
    end
    KPT_III = exp(-DG0_III/R0TP);
elseif phi>1  % case rich mixtures
    if (x == 0) && (y ~= 0) && M.L_minor>0 % if there are only hydrogens (H)
        for n = M.L_minor:-1:1
            DG0_VI(n) = (species_g0_new(minor_products{n},TP,strThProp) ...
                -C.alpha(n) * g_CO2 ...
                -(C.gamma(n)-2*C.alpha(n)) * g_H2O)*1000;
        end
        KPT_VI = exp(-DG0_VI/R0TP);
        DNfactor_VI = 1-C.alpha-C.beta/2-C.omega/2;
    elseif ((x ~= 0) && (y == 0) && M.L_minor>0 && phi < phi_c) && ~FLAG_SOOT% if there are only carbons (C)
        for n = M.L_minor:-1:1
            DG0_V(n) = (species_g0_new(minor_products{n},TP,strThProp) ...
                -(C.gamma(n)-C.alpha(n)-C.beta(n)/2) * g_CO2 ...
                -(C.beta(n)/2) * g_H2O ...
                -(2*C.alpha(n)-C.gamma(n)+C.beta(n)/2) * g_CO)*1000;
        end
        KPT_V = exp(-DG0_V/R0TP);
        DNfactor_V = 1-C.alpha-C.beta/2-C.omega/2;
    elseif phi < phi_c*TN.factor_c && ~FLAG_SOOT% general case of rich mixtures with hydrogens (H) and carbons (C)
        for n = M.L_minor:-1:1
            DG0_V(n) = (species_g0_new(minor_products{n},TP,strThProp) ...
                -(C.gamma(n)-C.alpha(n)-C.beta(n)/2) * g_CO2 ...
                -(C.beta(n)/2) * g_H2O ...
                -(2*C.alpha(n)-C.gamma(n)+C.beta(n)/2) * g_CO)*1000;
        end
        if M.L_minor>0
            KPT_V = exp(-DG0_V/R0TP);
            DNfactor_V = 1-C.alpha-C.beta/2-C.omega/2;
        end
        DG0_IV = (g_CO+g_H2O-g_CO2)*1000;
        KPT_IV = exp(-DG0_IV/R0TP);
    elseif phi >= phi_c*TN.factor_c || FLAG_SOOT
        for n = M.L_minor:-1:1
            DG0_V(n) = (species_g0_new(minor_products{n},TP,strThProp) ...
                -(C.gamma(n)-C.alpha(n)-C.beta(n)/2) * g_CO2 ...
                -(C.beta(n)/2) * g_H2O ...
                -(2*C.alpha(n)-C.gamma(n)+C.beta(n)/2) * g_CO)*1000;
        end
        if M.L_minor>0
            KPT_V = exp(-DG0_V/R0TP);
            DNfactor_V = 1-C.alpha-C.beta/2-C.omega/2;
        end
        
        DG0_IV = (g_CO+g_H2O-g_CO2)*1000;
        KPT_IV = exp(-DG0_IV/R0TP);
        DG0_VII = (g_CO2-2*g_CO)*1000;
        KPT_VII = exp(-DG0_VII/(C.R0*TP));
        mu = KPT_VII/DNfactor;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHEMICAL EQUILIBRIUM CALCULATIONS
if phi <= 1 % case of lean-to-stoichiometric mixtures
    [N_IC,DeltaNP] = incomplete_phi_1(N_IC,E,S,M,C);  
elseif phi > 1 && (x == 0) && (y ~= 0) %  case of rich mixtures if there are only hydrogens (H)
    [N_IC,DeltaNP] = incomplete_phi_2(N_IC,E,S,M,C);   
elseif ((x ~= 0) && (y == 0) && phi < phi_c) && ~FLAG_SOOT % if there are only carbons (C)
    [N_IC,DeltaNP] = incomplete_phi_3(N_IC,E,S,M,C);     
elseif phi < phi_c*TN.factor_c && ~FLAG_SOOT % general case of rich mixtures with hydrogens (H) and carbons (C) and without soot
    [N_IC,DeltaNP] = incomplete_phi_4(N_IC,E,S,M,C); 
elseif phi >= phi_c*TN.factor_c || FLAG_SOOT% rich mixtures with soot
    if (x ~= 0) && (y == 0) % if there are only carbons (C)
        [N_IC,DeltaNP] = incomplete_phi_5(N_IC,E,S,M,C);
    else % general case
        [N_IC,DeltaNP] = incomplete_phi_6(N_IC,E,S,M,C);
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
    function [N_IC,DeltaNP] = incomplete_phi_1(N_IC,E,S,M,C)
    DeltaNP = 1;
    while abs(DeltaNP/NP) > C.tolN
        % Initialization of the product matrix for incomplete combustion (IC)
        NP_old   = NP;
        N_IC_old = N_IC;
        N_IC = M0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % In lean combustion the product mixture always contains O2, CO2,
        % H2O, and N2 (in fuel-air combustion). The number of moles of
        % CO and H2 can be calculated from the equilibrium condition for
        % the reactions
        %
        %                  CO2 <-I -> CO+(1/2) O2              [KPT_I]
        %                  H2O <-II-> H2+(1/2) O2             [KPT_II]
        %
        % For the remaining minor species, we calculate the number of
        % moles from the equilibrium condition for the generalized
        % reaction
        %
        % C.alpha CO2+(C.beta/2) H2O+(C.gamma/2-C.alpha-C.beta/4) O2
        %+(C.omega/2) N2 <-III-> C_alpha H_beta O_gamma N_omega [KPT_III]
        %
        % Note that reactions*i*and II are particular forms of the more
        % general reaction III
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        NCO    = NCO2/NO2^(1/2)*DNfactor^(1/2)*KPT_I;
        NCO_old = N_IC_old(S.idx_CO,1);
        
        if NCO_old ~=0
            NCO = NCO_old+relax*(NCO-NCO_old);
        end
        N_IC(S.idx_CO,1) = NCO;
        
        NH2     = NH2O/NO2^(1/2)*DNfactor^(1/2)*KPT_II;
        NH2_old = N_IC_old(S.idx_H2,1);
        if NH2_old ~=0
            NH2 = NH2_old+relax*(NH2-NH2_old);
        end
        N_IC(S.idx_H2,1) = NH2;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Determination of the number of moles of the minor species
        % from the equilibrium condition for the above reaction
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if M.L_minor>0 
            Ni     = KPT_III.* NCO2.^C.alpha.* NH2O.^(C.beta/2).* NN2.^(C.omega/2) ...
                .* NO2.^(C.gamma/2-C.alpha-C.beta/4).* DNfactor.^DNfactor_III;
%           Ni = exp(log(KPT_III)+ log(NCO2).*C.alpha+ log(NH2O).*(C.beta/2)+ log(NN2).*(C.omega/2) ...
%                +log(NO2).*(C.gamma/2-C.alpha-C.beta/4)+ log(DNfactor).*DNfactor_III); 
%             for n =M.L_minor:-1:1
%                 Ni_old = N_IC_old(M.idx_minor(n),1);
%                 if Ni_old ~=0
%                     Ni(n) = Ni_old+relax*(Ni(n)-Ni_old);
%                 end
%                 N_IC(M.idx_minor(n),[strThProp.(minor_products{n}).Element_matrix(1,:),1]) = [Ni(n)*strThProp.(minor_products{n}).Element_matrix(2,:),Ni(n)];
%             end
            Ni_old = N_IC_old(M.idx_minor,1)';
            aux = find(Ni_old~=0);
            if ~isempty(aux)
                Ni(aux) = Ni_old(aux)+relax*(Ni(aux)-Ni_old(aux));
            end
            N_IC(M.idx_minor,1) = Ni;
        end
        % Correction of the number of moles of CO2, H2O, O2 and N2 from atom
        % conservation equations
        NCO2_old = NCO2;
        NCO2 = NCO2_0-sum(N_IC(:,1).*A0(:,E.ind_C)); % C-atom conservation
        if TP > 3000, NCO2 = NCO2_old+relax*(NCO2-NCO2_old); end
        if NCO2 < 0, NCO2 = 0.75*NCO2_old; end
        
        NH2O_old = NH2O;
        NH2O = NH2O_0-sum(N_IC(:,1).*A0(:,E.ind_H))/2; % H-atom conservation
        if TP > 3000, NH2O = NH2O_old+relax*(NH2O-NH2O_old); end
        if NH2O < 0, NH2O = 0.75*NH2O_old; end
        
        N_IC(S.idx_CO2,1)= NCO2;
        N_IC(S.idx_H2O,1)= NH2O;
        
        NO2_old = NO2;
        NO2  = NO2_0+NCO2_0+NCO_0/2+NH2O_0/2-sum(N_IC(:,1).*A0(:,E.ind_O))/2; % O-atom conservation
        if TP > 3000, NO2 = NO2_old+relax*(NO2-NO2_old); end
        if NO2 < 0, NO2 = 0.75*NO2_old; end
        
        N_IC(S.idx_O2,1) = NO2;
        NN2_old = NN2;
        NN2  = NN2_0 -sum(N_IC(:,1).*A0(:,E.ind_N))/2; % N-atom conservation
        if TP > 3000, NN2 = NN2_old+relax*(NN2-NN2_old); end
        if NN2 < 0, NN2 = 0.75*NN2_old; end
        
        N_IC(S.idx_N2,1) = NN2;
        N_IC(S.idx_He,1) = NHe;
        N_IC(S.idx_Ar,1) = NAr;
        % Overall number of moles in the product species
        NP = sum(N_IC(:,1).*(1-N_IC(:,10)));
%         NP = sum(N_IC(:,1));
        DeltaNP = norm([NP-NP_old, x-sum(N_IC(:,1).*A0(:,E.ind_C)), y-sum(N_IC(:,1).*A0(:,E.ind_H)), z-sum(N_IC(:,1).*A0(:,E.ind_O)), w-sum(N_IC(:,1).*A0(:,E.ind_N))]);
        % relax = abs(relax*(1-DeltaNP/NP));
    end
    end % PHI <= 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [N_IC,DeltaNP] = incomplete_phi_2(N_IC,E,S,M,C)
        DeltaNP = 1;
        while abs(DeltaNP/NP) > C.tolN
            % Initialization of the product matrix for incomplete combustion (IC)
            NP_old   = NP;
            N_IC_old = N_IC;
            
            N_IC = M0;
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
            % +(C.omega/2) N2 <-VI-> C_alpha H_beta O_gamma N_omega [KPT_VI]
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
            if M.L_minor>0
%                 Ni     = KPT_VI.* NCO2.^C.alpha.* NH2O.^(C.gamma-2*C.alpha) ...
%                     .* NH2.^(C.beta/2-C.gamma(n)+2*C.alpha) .* NN2.^(C.omega/2) ...
%                     .* DNfactor.^DNfactor_VI;
                Ni     = exp(log(KPT_VI)+log(NCO2).*(C.alpha)+log(NH2O).*(C.gamma-2*C.alpha) ...
                    +log(NH2).*(C.beta/2-C.gamma(n)+2*C.alpha)+log(NN2).*(C.omega/2) ...
                    +log(DNfactor).*DNfactor_VI);
                Ni_old = N_IC_old(M.idx_minor,1)';
                aux = find(Ni_old~=0);
                if ~isempty(aux)
                    Ni(aux) = Ni_old(aux)+relax*(Ni(aux)-Ni_old(aux));
                end
                N_IC(M.idx_minor,1) = Ni;
            end
            NO2     = DNfactor*(KPT_II*NH2O/NH2)^2;
            NO2_old = N_IC_old(S.idx_O2,1);
            if NO2_old ~=0
                NO2 = NO2_old+relax*(NO2-NO2_old);
            end
            
            % Correction of the number of moles of H2, H2O and N2 from
            % atom conservation equations
            
            NH2O = NH2O_0-2*(NO2+sum(N_IC(:,1).*A0(:,E.ind_O))/2);
            NH2  = NH2_0 +2*(NO2+sum(N_IC(:,1).*A0(:,E.ind_O))/2-sum(N_IC(:,1).*A0(:,E.ind_H))/4);
            NN2  = NN2_0 -sum(N_IC(:,1).*A0(:,E.ind_N))/2; % N-atom conservation

            N_IC(S.idx_H2O,1)= NH2O;
            N_IC(S.idx_H2,1) = NH2;
            N_IC(S.idx_O2,1) = NO2;
            N_IC(S.idx_N2,1) = NN2;
            N_IC(S.idx_He,1) = NHe;
            N_IC(S.idx_Ar,1) = NAr;
            % Overall number of moles in the product species
            NP = sum(N_IC(:,1).*(1-N_IC(:,10)));
%             NP = sum(N_IC(:,1));
            DeltaNP = NP-NP_old;
            % relax = abs(relax*(1-DeltaNP/NP));
        end
    end % 1 < PHI < PHI_C AND ~H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [N_IC,DeltaNP] = incomplete_phi_3(N_IC,E,S,M,C)
        DeltaNP = 1;
        while abs(DeltaNP/NP) > C.tolN
            % Initialization of the product matrix for incomplete combustion (IC)
            NP_old   = NP;
            N_IC_old = N_IC;
            
            N_IC = M0;
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
            %                 <-V-> C_alpha H_beta O_gamma N_omega   [KPT_V]
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
            if M.L_minor>0
%                 Ni     = KPT_V .* NCO2.^(C.gamma-C.alpha-C.beta/2) .* NH2O.^(C.beta/2) ...
%                     .* NCO.^(C.beta/2-C.gamma+2*C.alpha) .* NN2.^(C.omega/2) ...
%                     .* DNfactor.^DNfactor_V;
                Ni    = exp(log(KPT_V)+(C.gamma-C.alpha-C.beta/2).*log(NCO2)+(C.beta/2).*log(NH2O) ...
                       +log(NCO).*(2*C.alpha-C.gamma+C.beta/2)+log(NN2).*(C.omega/2) ...
                       +log(DNfactor).*DNfactor_V);
                Ni_old = N_IC_old(M.idx_minor,1)';
                aux = find(Ni_old~=0);
                if ~isempty(aux)
                    Ni(aux) = Ni_old(aux)+relax*(Ni(aux)-Ni_old(aux));
                end
                N_IC(M.idx_minor,1) = Ni;
            end
            NO2     = DNfactor*(KPT_I*NCO2/NCO)^2;
            NO2_old = N_IC_old(S.idx_O2,1);
            if NO2_old ~=0
                NO2 = NO2_old+relax*(NO2-NO2_old);
            end
            
            % Correction of the number of moles of H2, H2O and N2 from
            % atom conservation equations
            
            NCO2 = NCO2_0-2*(NO2-3*sum(N_IC(:,1).*A0(:,E.ind_C))/2+sum(N_IC(:,1).*A0(:,E.ind_O))/2);
            NCO  = NCO_0 +2*(NO2-sum(N_IC(:,1).*A0(:,E.ind_C))+sum(N_IC(:,1).*A0(:,E.ind_O))/2);
            NN2  = NN2_0 -sum(N_IC(:,1).*A0(:,E.ind_N))/2; % N-atom conservation
            
            N_IC(S.idx_CO,1) = NCO;
            N_IC(S.idx_CO2,1)= NCO2;
            N_IC(S.idx_O2,1) = NO2;
            N_IC(S.idx_N2,1) = NN2;
            N_IC(S.idx_He,1) = NHe;
            N_IC(S.idx_Ar,1) = NAr;
            % Overall number of moles in the product species
            NP = sum(N_IC(:,1).*(1-N_IC(:,10)));
%             NP = sum(N_IC(:,1));
            DeltaNP = NP-NP_old;
            % relax = abs(relax*(1-DeltaNP/NP));
        end
    end % 1 < PHI < PHI_C AND ~C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [N_IC,DeltaNP] = incomplete_phi_4(N_IC,E,S,M,C)
                DeltaNP = 1;
        while abs(DeltaNP/NP) > C.tolN && it<itMax 
            it = it+1;
            % Initialization of the product matrix for incomplete combustion (IC)
            NP_old   = NP;
            N_IC_old = N_IC;
            
            N_IC = M0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % In rich combustion the product mixture always contains CO2, CO,
            % H2O, H2 and N2 (in fuel-air combustion). The number of moles of
            % these species must be calculated combining the 4 atom
            % conservation equations with the equilibrium condition for
            % the inverse water-gas shift reaction
            %
            %                  CO2+H2 <-IV-> H2O+CO             [KPT_IV]
            %
            % For the remaining minor species, we calculate the number of
            % moles from the equilibrium condition for the generalized
            % reactions
            %
            % (C.gamma-C.alpha-C.beta/2) CO2+(C.beta/2) H2O
            %    +(2*C.alpha-C.gamma+C.beta/2) CO+(C.omega/2) N2
            %                 <-V-> C_alpha H_beta O_gamma N_omega   [KPT_V]
            %
            % or
            %
            % C.alpha CO2+(C.gamma-2*C.alpha) H2O+(C.beta/2-C.gamma+2*C.alpha) H2
            % +(C.omega/2) N2 <-VI-> C_alpha H_beta O_gamma N_omega [KPT_VI]
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
            % low Hydrogen:    CO2 <-I -> CO+(1/2) O2              [KPT_I]
            % low Carbon:      H2O <-II-> H2+(1/2) O2             [KPT_II]
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Determination of the number of moles of the minor species from
            % the equilibrium condition for the above reaction
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if M.L_minor>0 
%                 Ni     = KPT_V .* NCO2.^(C.gamma-C.alpha-C.beta/2) .* NH2O.^(C.beta/2) ...
%                     .* NCO.^(2*C.alpha-C.gamma+C.beta/2) .* NN2.^(C.omega/2) ...
%                     .* DNfactor.^DNfactor_V;
                Ni    = exp(log(KPT_V)+(C.gamma-C.alpha-C.beta/2).*log(NCO2)+(C.beta/2).*log(NH2O) ...
                       +log(NCO).*(2*C.alpha-C.gamma+C.beta/2)+log(NN2).*(C.omega/2) ...
                       +log(DNfactor).*DNfactor_V);
                Ni_old = N_IC_old(M.idx_minor,1)';
                aux = find(Ni_old~=0);
                if ~isempty(aux)
                    Ni(aux) = Ni_old(aux)+relax*(Ni(aux)-Ni_old(aux));
                end
                N_IC(M.idx_minor,1) = Ni;
            end
            
            NO2     = DNfactor*(KPT_I*NCO2/NCO)^2;
            NO2_old = N_IC_old(S.idx_O2,1);
            if NO2_old ~=0
                NO2 = NO2_old+relax*(NO2-NO2_old);
            end
            
            % Correction of the number of moles of CO, H2, CO2, H2O and N2 from
            % atom conservation equations
            
            a = NCO2_0+NCO_0-sum(N_IC(:,1).*A0(:,E.ind_C));
            b = NH2O_0+NH2_0-sum(N_IC(:,1).*A0(:,E.ind_H))/2;
            c = NCO_0+NH2_0+2*(NO2-sum(N_IC(:,1).*A0(:,E.ind_C))-sum(N_IC(:,1).*A0(:,E.ind_H))/4+sum(N_IC(:,1).*A0(:,E.ind_O))/2);
            d = b-c;
            
            NCO2_old = NCO2;
            NCO_old  = NCO;
            NH2O_old = NH2O;
            NH2_old  = NH2;
            NN2_old  = NN2;
            
            NCO  = (1/2)*((a+c)*KPT_IV+d-sqrt((a-c)^2*KPT_IV^2+(2*(2*a*c+d*(a+c)))*KPT_IV+d^2))/(KPT_IV-1);
            NH2  = c-NCO;
            NH2O = d+NCO;
            NCO2 = a-NCO;
            NN2  = NN2_0 -sum(N_IC(:,1).*A0(:,E.ind_N))/2; % N-atom conservation
            
            if TP > 3000
                NCO2 = NCO2_old+0.1*relax*(NCO2-NCO2_old);
                NCO  = NCO_old +0.1*relax*(NCO -NCO_old);
                NH2O = NH2O_old+0.1*relax*(NH2O-NH2O_old);
                NH2  = NH2_old +0.1*relax*(NH2 -NH2_old);
                NN2  = NN2_old +0.1*relax*(NN2 -NN2_old);
            end
             
            if NCO2 < 0, NCO2 = 0.75*NCO2_old; end
            if NCO  < 0, NCO  = 0.75*NCO_old ; end
            if NH2O < 0, NH2O = 0.75*NH2O_old; end
            if NH2  < 0, NH2  = 0.75*NH2_old ; end
            if NN2  < 0, NN2  = 0.75*NN2_old ; end
            
            N_IC(S.idx_CO,1) = NCO;
            N_IC(S.idx_CO2,1)= NCO2;
            N_IC(S.idx_H2O,1)= NH2O;
            N_IC(S.idx_H2,1) = NH2;
            N_IC(S.idx_O2,1) = NO2;
            N_IC(S.idx_N2,1) = NN2;
            N_IC(S.idx_He,1) = NHe;
            N_IC(S.idx_Ar,1) = NAr;
            
            % Overall number of moles in the product species
            
            NP = sum(N_IC(:,1).*(1-N_IC(:,10)));
%             NP = sum(N_IC(:,1));
            DeltaNP = norm([NP-NP_old, x-sum(N_IC(:,1).*A0(:,E.ind_C)), y-sum(N_IC(:,1).*A0(:,E.ind_H)), z-sum(N_IC(:,1).*A0(:,E.ind_O)), w-sum(N_IC(:,1).*A0(:,E.ind_N))]);
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
    function [N_IC,DeltaNP] = incomplete_phi_5(N_IC,E,S,M,C)
        error('NOT IMPLEMENTED YET');
%         DeltaNP = 1;
%         while abs(DeltaNP/NP) > C.tolN
%             % Initialization of the product matrix for incomplete combustion (IC)
%             NP_old   = NP;
%             N_IC_old = N_IC;
%             
%             N_IC = M0;
%         end
    end % PHI <= 1 NOT IMPLEMENTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [N_IC,DeltaNP] = incomplete_phi_6(N_IC,E,S,M,C)
        DeltaNP = 1;
        while abs(DeltaNP/NP) > C.tolN && it<itMax 
            it = it+1;
            % Initialization of the product matrix for incomplete combustion (IC)
            NP_old   = NP;
            N_IC_old = N_IC;
                      
            N_IC = M0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % In rich combustion the product mixture always contains CO2, CO,
            % H2O, H2 and N2 (in fuel-air combustion). The number of moles of
            % these species must be calculated combining the 4 atom
            % conservation equations with the equilibrium condition for
            % the inverse water-gas shift reaction
            %
            %                  CO2+H2 <-IV-> H2O+CO             [KPT_IV]
            %
            % For the remaining minor species, we calculate the number of
            % moles from the equilibrium condition for the generalized
            % reactions
            %
            % (C.gamma-C.alpha-C.beta/2) CO2+(C.beta/2) H2O
            %    +(2*C.alpha-C.gamma+C.beta/2) CO+(C.omega/2) N2
            %                 <-V-> C_alpha H_beta O_gamma N_omega   [KPT_V]
            %
            % or
            %
            % C.alpha CO2+(C.gamma-2*C.alpha) H2O+(C.beta/2-C.gamma+2*C.alpha) H2
            % +(C.omega/2) N2 <-VI-> C_alpha H_beta O_gamma N_omega [KPT_VI]
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
            % low Hydrogen:    CO2 <-I -> CO+(1/2) O2              [KPT_I]
            % low Carbon:      H2O <-II-> H2+(1/2) O2             [KPT_II]
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
%                 if M.L_minor>0
%                     Ni    = KPT_VIII .* NCO.^(C.gamma) .* NCgr.^(C.alpha-C.gamma) ...
%                         .* NH2.^(C.beta/2) .* NN2.^(C.omega/2) ...
%                         .* DNfactor.^DNfactor_VIII;
%                     for n = M.L_minor:-1:1
%                         Ni_old = N_IC_old(M.idx_minor(n),1);
%                         if Ni_old ~=0
%                             Ni(n) = Ni_old+relax*(Ni(n)-Ni_old);
%                         end
%                         N_IC(M.idx_minor(n),[strThProp.(minor_products{n}).Element_matrix(1,:),1]) = [Ni(n)*strThProp.(minor_products{n}).Element_matrix(2,:),Ni(n)];
%                     end
%                 end
%             else
                if M.L_minor>0
                    Ni    = KPT_V .* NCO2.^(C.gamma-C.alpha-C.beta/2) .* NH2O.^(C.beta/2) ...
                        .* NCO.^(2*C.alpha-C.gamma+C.beta/2) .* NN2.^(C.omega/2) ...
                        .* DNfactor.^DNfactor_V;
%                     Ni    = exp(log(KPT_V)+(C.gamma-C.alpha-C.beta/2).*log(NCO2)+(C.beta/2).*log(NH2O) ...
%                        +log(NCO).*(2*C.alpha-C.gamma+C.beta/2)+log(NN2).*(C.omega/2) ...
%                        +log(DNfactor).*DNfactor_V);

                    Ni_old = N_IC_old(M.idx_minor,1)';
                    aux = find(Ni_old~=0);
                    if ~isempty(aux)
                        Ni(aux) = Ni_old(aux)+relax*(Ni(aux)-Ni_old(aux));
                    end
                    N_IC(M.idx_minor,1) = Ni;
                end
%             end
%             NO2     = mean([DNfactor*(KPT_I*NCO2/NCO)^2, DNfactor*(KPT_II*NH2O/NH2)^2]);
            NO2     = DNfactor*(KPT_I*NCO2/NCO)^2;
%             NO2 = 0;
            NO2_old = N_IC_old(S.idx_O2,1);
            if NO2_old ~= 0
                NO2 = NO2_old+relax*(NO2-NO2_old);
            end
            
            % Correction of the number of moles of CO, H2, CO2, H2O and N2 from
            % atom conservation equations
            
            a = NCO2_0+NCO_0+NCgr_0-sum(N_IC(:,1).*A0(:,E.ind_C));
            b = NH2O_0+NH2_0-sum(N_IC(:,1).*A0(:,E.ind_H))/2;
            c = 2*NCO2_0+NCO_0+NH2O_0-2*NO2-sum(N_IC(:,1).*A0(:,E.ind_O));
            
            %%% NEW
%             a0 = 4*a^3*KPT_VII^2+4*a^2*b*KPT_VII^2+a*b^2*KPT_VII^2-4*a^2*c*KPT_VII^2-2*a*b*c*KPT_VII^2+a*c^2*KPT_VII^2+(a*KPT_IV)/DNfactor^2-(c*KPT_IV)/DNfactor^2-(a*KPT_IV^2)/DNfactor^2+(c*KPT_IV^2)/DNfactor^2+(a^2*KPT_VII)/DNfactor+(a*b*KPT_VII)/DNfactor-(a*c*KPT_VII)/DNfactor+(4*a^2*KPT_IV*KPT_VII)/DNfactor-(4*a*c*KPT_IV*KPT_VII)/DNfactor-(b*c*KPT_IV*KPT_VII)/DNfactor+(c^2*KPT_IV*KPT_VII)/DNfactor-(4*a^2*KPT_IV^2*KPT_VII)/DNfactor+(4*a*c*KPT_IV^2*KPT_VII)/DNfactor-(c^2*KPT_IV^2*KPT_VII)/DNfactor;
%             a1 = 12*a^2*KPT_VII^2+8*a*b*KPT_VII^2+b^2*KPT_VII^2-8*a*c*KPT_VII^2-2*b*c*KPT_VII^2+c^2*KPT_VII^2+KPT_IV/DNfactor^2-KPT_IV^2/DNfactor^2+(2*a*KPT_VII)/DNfactor+(b*KPT_VII)/DNfactor-(c*KPT_VII)/DNfactor+(8*a*KPT_IV*KPT_VII)/DNfactor-(4*c*KPT_IV*KPT_VII)/DNfactor-(8*a*KPT_IV^2*KPT_VII)/DNfactor+(4*c*KPT_IV^2*KPT_VII)/DNfactor;
%             a2 = (-12*a*KPT_VII-4*b*KPT_VII+4*c*KPT_VII-1/DNfactor-(4*KPT_IV)/DNfactor+(4*KPT_IV^2)/DNfactor)*KPT_VII;
%             a3 = 4*KPT_VII^2;
            %%% NEW
            
            NCO2_old = NCO2;
            NCO_old  = NCO;
            NCgr_old = NCgr;
            NH2O_old = NH2O;
            NH2_old  = NH2;
            NN2_old  = NN2;
                       
            NCgr = real(1/24*(24*a+(4*(2+KPT_IV))/(KPT_IV*mu)+(2*(1+1i*sqrt(3))*(-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu)))/(KPT_IV*(8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))+(2*1i*(1i+sqrt(3))*(8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))/(KPT_IV*mu^2)-(1/(6*KPT_IV^2*mu^3))*((2*(2+KPT_IV)*mu+((1+1i*sqrt(3))*mu^2*(-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu)))/(8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3)+1i*(1i+sqrt(3))*(8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))^2)));  
%             NCgr = real(1/24*(24*a+(4*(2+KPT_IV))/(KPT_IV*mu)+(2*(1-1i*sqrt(3))*(-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu)))/(KPT_IV*(8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))-(2*1i*(-1i+sqrt(3))*(8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))/(KPT_IV*mu^2)-(1/(6*KPT_IV^2*mu^3))*((-2*(2+KPT_IV)*mu+(1i*(1i+sqrt(3))*mu^2*(-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu)))/(8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3)+(1+1i*sqrt(3))*(8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))^2)));  
            %%% NEW
%             NCgr = real(-(a2/(3*a3))+(2^(1/3)*(-a2^2+3*a1*a3))/(3*a3*(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3))-(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3)/(3*2^(1/3)*a3));
%             NCgr = real(-(a2/(3*a3))-((1-1i*sqrt(3))*(-a2^2+3*a1*a3))/(3*2^(2/3)*a3*(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3))+((1+1i*sqrt(3))*(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3))/(6*2^(1/3)*a3));
%             NCgr = real(-(a2/(3*a3))-((1+1i*sqrt(3))*(-a2^2+3*a1*a3))/(3*2^(2/3)*a3*(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3))+((1-1i*sqrt(3))*(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3))/(6*2^(1/3)*a3));
            %%% NEW
            NCO2 = (1+2*a*mu-2*mu*NCgr-sqrt(1+4*a*mu-4*mu*NCgr))/(2*mu);
            %%% NEW 2 
%             a0 = -2*KPT_IV/DNfactor-KPT_VII*b+2*KPT_VII*c;
%             a1 = 2*KPT_VII+4*KPT_IV*KPT_VII;
%             a2 = 4*KPT_VII^2*DNfactor;
%             a3 = 2*KPT_IV*c/DNfactor;
%             NCO = real(-(a1/(3*a2))-(2^(1/3)*(-a1^2-3*a0*a2))/(3*a2*(-2*a1^3-9*a0*a1*a2+27*a2^2*a3+sqrt(-4*(a1^2+3*a0*a2)^3+(2*a1^3+9*a0*a1*a2-27*a2^2*a3)^2))^(1/3))+(-2*a1^3-9*a0*a1*a2+27*a2^2*a3+sqrt(-4*(a1^2+3*a0*a2)^3+(2*a1^3+9*a0*a1*a2-27*a2^2*a3)^2))^(1/3)/(3*2^(1/3)*a2));
%             NCO2 = NCO^2*mu;
            %%% NEW 2
            if ~isreal(NCO2)  
            	NCgr = real(1/6*(6*a+(2+KPT_IV)/(KPT_IV*mu)+(4-2*KPT_IV+KPT_IV^2*(1-6*b*mu+6*c*mu))/(KPT_IV*(8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))+(1/(KPT_IV*mu^2))*((8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))-(1/(6*KPT_IV^2*mu^3))*((2*mu+KPT_IV*mu-(mu^2*(-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu)))/(8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3)+(8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))^2)));
                %%% NEW
%                 NCgr = real(-(a2/(3*a3))+(2^(1/3)*(-a2^2+3*a1*a3))/(3*a3*(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3))-(2*a2^3-9*a1*a2*a3-27*a0*a3^2+sqrt(-4*(a2^2-3*a1*a3)^3+(-2*a2^3+9*a1*a2*a3+27*a0*a3^2)^2))^(1/3)/(3*2^(1/3)*a3));
                %%% NEW
                NCO2 = (1+2*a*mu-2*mu*NCgr-sqrt(1+4*a*mu-4*mu*NCgr))/(2*mu);
            end
            NCO = (-1+sqrt(1+4*a*mu-4*mu*NCgr))/(2*mu);
            NH2O = c-2*NCO2-NCO;
            NH2 = b-NH2O;
            %%%%%%%%%%%%%%%%%%%%%%%%
%             NH2O = real((1/(6*KPT_IV))*(2*(KPT_IV*(2*b+c-2*mu2)+mu2)+(2*2^(1/3)*(b^2*KPT_IV^2+(c*KPT_IV+mu2-2*KPT_IV*mu2)^2+b*KPT_IV*(mu2-2*KPT_IV*(c+4*mu2))))/(-2*b^3*KPT_IV^3+3*b^2*KPT_IV^2*(2*KPT_IV*(c-5*mu2)-mu2)+2*(c*KPT_IV+mu2-2*KPT_IV*mu2)^3-3*b*KPT_IV*(c*KPT_IV+mu2-2*KPT_IV*mu2)*(2*c*KPT_IV+(-1+8*KPT_IV)*mu2)+3*sqrt(3)*sqrt(-b^2*KPT_IV^2*mu2*(-8*b^3*KPT_IV^4+(8*c*KPT_IV^2+mu2)*(c*KPT_IV+mu2-2*KPT_IV*mu2)^2+2*b*KPT_IV*(-12*c^2*KPT_IV^3+c*KPT_IV*(-1+2*KPT_IV-40*KPT_IV^2)*mu2+(1+4*(-2+KPT_IV)*KPT_IV)*mu2^2)+b^2*KPT_IV^2*(mu2+4*KPT_IV*(6*c*KPT_IV+(-5+KPT_IV)*mu2)))))^(1/3)+2^(2/3)*(-2*b^3*KPT_IV^3+3*b^2*KPT_IV^2*(2*KPT_IV*(c-5*mu2)-mu2)+2*(c*KPT_IV+mu2-2*KPT_IV*mu2)^3-3*b*KPT_IV*(c*KPT_IV+mu2-2*KPT_IV*mu2)*(2*c*KPT_IV+(-1+8*KPT_IV)*mu2)+3*sqrt(3)*sqrt(-b^2*KPT_IV^2*mu2*(-8*b^3*KPT_IV^4+(8*c*KPT_IV^2+mu2)*(c*KPT_IV+mu2-2*KPT_IV*mu2)^2+2*b*KPT_IV*(-12*c^2*KPT_IV^3+c*KPT_IV*(-1+2*KPT_IV-40*KPT_IV^2)*mu2+(1+4*(-2+KPT_IV)*KPT_IV)*mu2^2)+b^2*KPT_IV^2*(mu2+4*KPT_IV*(6*c*KPT_IV+(-5+KPT_IV)*mu2)))))^(1/3)));
% %             NH2O = real((1/(12*KPT_IV))*(4*(KPT_IV*(2*b+c-2*mu2)+mu2)-(2*(-1i)*2^(1/3)*(-(-1i)+sqrt(3))*(b^2*KPT_IV^2+(c*KPT_IV+mu2-2*KPT_IV*mu2)^2+b*KPT_IV*(mu2-2*KPT_IV*(c+4*mu2))))/(-2*b^3*KPT_IV^3+3*b^2*KPT_IV^2*(2*KPT_IV*(c-5*mu2)-mu2)+2*(c*KPT_IV+mu2-2*KPT_IV*mu2)^3-3*b*KPT_IV*(c*KPT_IV+mu2-2*KPT_IV*mu2)*(2*c*KPT_IV+(-1+8*KPT_IV)*mu2)+3*sqrt(3)*sqrt(-b^2*KPT_IV^2*mu2*(-8*b^3*KPT_IV^4+(8*c*KPT_IV^2+mu2)*(c*KPT_IV+mu2-2*KPT_IV*mu2)^2+2*b*KPT_IV*(-12*c^2*KPT_IV^3+c*KPT_IV*(-1+2*KPT_IV-40*KPT_IV^2)*mu2+(1+4*(-2+KPT_IV)*KPT_IV)*mu2^2)+b^2*KPT_IV^2*(mu2+4*KPT_IV*(6*c*KPT_IV+(-5+KPT_IV)*mu2)))))^(1/3)+(-1i)*2^(2/3)*((-1i)+sqrt(3))*(-2*b^3*KPT_IV^3+3*b^2*KPT_IV^2*(2*KPT_IV*(c-5*mu2)-mu2)+2*(c*KPT_IV+mu2-2*KPT_IV*mu2)^3-3*b*KPT_IV*(c*KPT_IV+mu2-2*KPT_IV*mu2)*(2*c*KPT_IV+(-1+8*KPT_IV)*mu2)+3*sqrt(3)*sqrt(-b^2*KPT_IV^2*mu2*(-8*b^3*KPT_IV^4+(8*c*KPT_IV^2+mu2)*(c*KPT_IV+mu2-2*KPT_IV*mu2)^2+2*b*KPT_IV*(-12*c^2*KPT_IV^3+c*KPT_IV*(-1+2*KPT_IV-40*KPT_IV^2)*mu2+(1+4*(-2+KPT_IV)*KPT_IV)*mu2^2)+b^2*KPT_IV^2*(mu2+4*KPT_IV*(6*c*KPT_IV+(-5+KPT_IV)*mu2)))))^(1/3)));
%             NH2 = b-NH2O;
%             if NH2 < 0 || ~isreal(NH2)
%                 NH2O = real((1/(12*KPT_IV))*(4*(KPT_IV*(2*b+c-2*mu2)+mu2)-(2*1i*2^(1/3)*(-1i+sqrt(3))*(b^2*KPT_IV^2+(c*KPT_IV+mu2-2*KPT_IV*mu2)^2+b*KPT_IV*(mu2-2*KPT_IV*(c+4*mu2))))/(-2*b^3*KPT_IV^3+3*b^2*KPT_IV^2*(2*KPT_IV*(c-5*mu2)-mu2)+2*(c*KPT_IV+mu2-2*KPT_IV*mu2)^3-3*b*KPT_IV*(c*KPT_IV+mu2-2*KPT_IV*mu2)*(2*c*KPT_IV+(-1+8*KPT_IV)*mu2)+3*sqrt(3)*sqrt(-b^2*KPT_IV^2*mu2*(-8*b^3*KPT_IV^4+(8*c*KPT_IV^2+mu2)*(c*KPT_IV+mu2-2*KPT_IV*mu2)^2+2*b*KPT_IV*(-12*c^2*KPT_IV^3+c*KPT_IV*(-1+2*KPT_IV-40*KPT_IV^2)*mu2+(1+4*(-2+KPT_IV)*KPT_IV)*mu2^2)+b^2*KPT_IV^2*(mu2+4*KPT_IV*(6*c*KPT_IV+(-5+KPT_IV)*mu2)))))^(1/3)+1i*2^(2/3)*(1i+sqrt(3))*(-2*b^3*KPT_IV^3+3*b^2*KPT_IV^2*(2*KPT_IV*(c-5*mu2)-mu2)+2*(c*KPT_IV+mu2-2*KPT_IV*mu2)^3-3*b*KPT_IV*(c*KPT_IV+mu2-2*KPT_IV*mu2)*(2*c*KPT_IV+(-1+8*KPT_IV)*mu2)+3*sqrt(3)*sqrt(-b^2*KPT_IV^2*mu2*(-8*b^3*KPT_IV^4+(8*c*KPT_IV^2+mu2)*(c*KPT_IV+mu2-2*KPT_IV*mu2)^2+2*b*KPT_IV*(-12*c^2*KPT_IV^3+c*KPT_IV*(-1+2*KPT_IV-40*KPT_IV^2)*mu2+(1+4*(-2+KPT_IV)*KPT_IV)*mu2^2)+b^2*KPT_IV^2*(mu2+4*KPT_IV*(6*c*KPT_IV+(-5+KPT_IV)*mu2)))))^(1/3)));
%                 NH2 = b-NH2O;
%             end
%             NCO2 = mu2*(NH2O/NH2)^2;
%             NCO = c-2*NCO2-NH2O;
%             NCgr = a-NCO2-NCO;
            %%%%%%%%%%%%%%%%%%%%%%%%
            NN2  = NN2_0-sum(N_IC(:,1).*A0(:,E.ind_N))/2; % N-atom conservation
            
            if TP > 3000
                NCO2 = NCO2_old+0.5*relax*(NCO2-NCO2_old);
                NCO  = NCO_old +0.5*relax*(NCO -NCO_old);
                NCgr = NCgr_old+0.5*relax*(NCgr-NCgr_old);
                NH2O = NH2O_old+0.5*relax*(NH2O-NH2O_old);
                NH2  = NH2_old +0.5*relax*(NH2 -NH2_old);
                NN2  = NN2_old +0.5*relax*(NN2 -NN2_old);
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
                    
                        N_minor_C_0 = sum(N_IC_old(:,1).*A0(:,E.ind_C))-NCO2_0-NCO_0-NCgr_0;
                        N_minor_H_0 = sum(N_IC_old(:,1).*A0(:,E.ind_H))-2*NH2O_0-2*NH2_0;
                        N_minor_O_0 = sum(N_IC_old(:,1).*A0(:,E.ind_O))-2*NO2_0-2*NCO2_0-NCO_0-NH2O_0;
                        N_minor_N_0 = sum(N_IC_old(:,1).*A0(:,E.ind_N))-2*NN2_0;
                        
                        N_minor_C_0_nswt = sum(N_IC_old(:,1).*A0(:,E.ind_C).*(1-N_IC(:,10)))-NCO2_0-NCO_0-NCgr_0;
                        N_minor_H_0_nswt = sum(N_IC_old(:,1).*A0(:,E.ind_H).*(1-N_IC(:,10)))-2*NH2O_0-2*NH2_0;
                        N_minor_O_0_nswt = sum(N_IC_old(:,1).*A0(:,E.ind_O).*(1-N_IC(:,10)))-2*NO2_0-2*NCO2_0-NCO_0-NH2O_0;
                        N_minor_N_0_nswt = sum(N_IC_old(:,1).*A0(:,E.ind_N).*(1-N_IC(:,10)))-2*NN2_0;
                    else
%                         NCO_0  = NCO;
%                         NCO2_0 = NCO2;
%                         NH2O_0 = NH2O;
%                         NH2_0  = NH2;
%                         NN2_0  = NN2;
%                         NO2_0  = NO2;
%                         NCgr_0 = NCgr;
                        
                        N_minor_C_0 = sum(N_IC(:,1).*A0(:,E.ind_C));
                        N_minor_H_0 = sum(N_IC(:,1).*A0(:,E.ind_H));
                        N_minor_O_0 = sum(N_IC(:,1).*A0(:,E.ind_O));
                        N_minor_N_0 = sum(N_IC(:,1).*A0(:,E.ind_N));
                        
                        N_minor_C_0_nswt = sum(N_IC(:,1).*A0(:,E.ind_C).*(1-N_IC(:,10)));
                        N_minor_H_0_nswt = sum(N_IC(:,1).*A0(:,E.ind_H).*(1-N_IC(:,10)));
                        N_minor_O_0_nswt = sum(N_IC(:,1).*A0(:,E.ind_O).*(1-N_IC(:,10)));
                        N_minor_N_0_nswt = sum(N_IC(:,1).*A0(:,E.ind_N).*(1-N_IC(:,10)));
                        
%                         NP_old = N_minor_C_0+N_minor_H_0+N_minor_O_0+N_minor_N_0...
%                             +NCO2_0+NCO_0+NH2O_0+NH2_0+NO2_0+NN2_0+NCgr_0...
%                             +NHe+NAr;
%                         NP_old = N_minor_C_0_nswt+N_minor_H_0_nswt+N_minor_O_0_nswt+N_minor_N_0_nswt...
%                             +NCO2_0+NCO_0+NH2O_0+NH2_0+NO2_0+NN2_0...
%                             +NHe+NAr;
                    end
                    
                    NP_old = N_minor_C_0_nswt+N_minor_H_0_nswt+N_minor_O_0_nswt+N_minor_N_0_nswt...
                            +y/4+0.5+(NH2_0+NCO_0+z+w)/2+NHe+NAr;
                    NP_old = 0.5*(NP+NP_old); 
                    
                    if strfind(PD.ProblemType,'P') == 2 % PD.ProblemType = 'TP', 'HP', or 'SP'
                        DNfactor = NP_old/pP;
                    elseif strfind(PD.ProblemType,'V') == 2 % PD.ProblemType = 'TV', 'EV', or 'SV'
                        DNfactor = (vP/1000)*1e5/R0TP; % vP is given in l
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
                if M.L_minor>0
%                     Ni    = KPT_V .* NCO2.^(C.gamma-C.alpha-C.beta/2) .* NH2O.^(C.beta/2) ...
%                         .* NCO.^(2*C.alpha-C.gamma+C.beta/2) .* NN2.^(C.omega/2) ...
%                         .* DNfactor.^DNfactor_V;
                    Ni    = exp(log(KPT_V)+(C.gamma-C.alpha-C.beta/2).*log(NCO2)+(C.beta/2).*log(NH2O) ...
                       +log(NCO).*(2*C.alpha-C.gamma+C.beta/2)+log(NN2).*(C.omega/2) ...
                       +log(DNfactor).*DNfactor_V);
                    Ni_old = N_IC_old(M.idx_minor,1)';
                    aux = find(Ni_old~=0);
                    if ~isempty(aux)
                        Ni(aux) = Ni_old(aux)+relax*(Ni(aux)-Ni_old(aux));
%                           Ni(aux) = exp(log(Ni_old(aux)) +relax*(log(Ni(aux)) -log(Ni_old(aux))));
                    end
                    N_IC(M.idx_minor,1) = Ni;
                end
                
                NO2     = DNfactor*(KPT_I*NCO2/NCO)^2;
                
                NO2_old = N_IC_old(S.idx_O2,1);
                if NO2_old ~= 0
                    NO2 = NO2_old+relax*(NO2-NO2_old);
%                     NO2 = exp(log(NO2_old) +relax*(log(NO2) -log(NO2_old)));
                end
                
                a = NCO2_0+NCO_0+NCgr_0+N_minor_C_0-sum(N_IC(:,1).*A0(:,E.ind_C));
                b = NH2O_0+NH2_0+N_minor_H_0/2-sum(N_IC(:,1).*A0(:,E.ind_H))/2;
                c = NCO_0+NH2_0+2*(NCgr_0+N_minor_C_0+N_minor_H_0/4-N_minor_O_0/2)+2*(NO2-NO2_0-sum(N_IC(:,1).*A0(:,E.ind_C))-sum(N_IC(:,1).*A0(:,E.ind_H))/4+sum(N_IC(:,1).*A0(:,E.ind_O))/2);

                d = b-c;
                
                NCO2_old = NCO2;
                NCO_old  = NCO;
                NH2O_old = NH2O;
                NH2_old  = NH2;
                NN2_old  = NN2;
                
                NCO  = (1/2)*((a+c)*KPT_IV+d-sqrt((a-c)^2*KPT_IV^2+(2*(2*a*c+d*(a+c)))*KPT_IV+d^2))/(KPT_IV-1);
                NCO2 = a-NCO;
                NH2  = c-NCO;
                NH2O = d+NCO;

                NN2  = NN2_0+N_minor_N_0/2-sum(N_IC(:,1).*A0(:,E.ind_N))/2; % N-atom conservation
                
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
                    NCO2 = NCO2_old+0.5*relax*(NCO2-NCO2_old);
                    NCO  = NCO_old +0.5*relax*(NCO -NCO_old);
                    NH2O = NH2O_old+0.5*relax*(NH2O-NH2O_old);
                    NH2  = NH2_old +0.5*relax*(NH2 -NH2_old);
                    NN2  = NN2_old +0.5*relax*(NN2 -NN2_old);
                end

                t = false;
            end
           
            
            if NCO2 < 0 || ~isreal(NCO2), NCO2 = 0.75*NCO2_old; end
            if NCO  < 0 || ~isreal(NCO), NCO  = 0.75*NCO_old ; end
            if NCgr < 0 || ~isreal(NCgr), NCgr = 0.01*NCgr_old ; end
            if NH2O < 0 || ~isreal(NH2O), NH2O = 0.75*NH2O_old; end
            if NH2  < 0 || ~isreal(NH2), NH2  = 0.75*NH2_old ; end
            if NN2  < 0 || ~isreal(NN2), NN2  = 0.75*NN2_old ; end
            
            N_IC(S.idx_Cgr,1)= NCgr;
            N_IC(S.idx_CO,1) = NCO;
            N_IC(S.idx_CO2,1)= NCO2;
            N_IC(S.idx_H2O,1)= NH2O;
            N_IC(S.idx_H2,1) = NH2;
            N_IC(S.idx_O2,1) = NO2;
            N_IC(S.idx_N2,1) = NN2;
            N_IC(S.idx_He,1) = NHe;
            N_IC(S.idx_Ar,1) = NAr;

            % Overall number of moles in the product species
            NP = sum(N_IC(:,1).*(1-N_IC(:,10)));
%             NP = sum(N_IC(:,1));
            DeltaNP = norm([NP-NP_old, x-sum(N_IC(:,1).*A0(:,E.ind_C)), y-sum(N_IC(:,1).*A0(:,E.ind_H)), z-sum(N_IC(:,1).*A0(:,E.ind_O)), w-sum(N_IC(:,1).*A0(:,E.ind_N))]);
            % relax = abs(relax*(1-DeltaNP/NP));
        end
    end % PHI >= PHI_C OR SOOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [N_IC,DeltaNP] = incomplete_phi_7(N_IC,E,S,M,C)
        DeltaNP = 1;
        while abs(DeltaNP/NP) > C.tolN && it<itMax 
            it = it+1;
            % Initialization of the product matrix for incomplete combustion (IC)
            NP_old   = NP;
            N_IC_old = N_IC;
                      
            N_IC = M0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % In rich combustion the product mixture always contains CO2, CO,
            % H2O, H2 and N2 (in fuel-air combustion). The number of moles of
            % these species must be calculated combining the 4 atom
            % conservation equations with the equilibrium condition for
            % the inverse water-gas shift reaction
            %
            %                  CO2+H2 <-IV-> H2O+CO             [KPT_IV]
            %
            % For the remaining minor species, we calculate the number of
            % moles from the equilibrium condition for the generalized
            % reactions
            %
            % (C.gamma-C.alpha-C.beta/2) CO2+(C.beta/2) H2O
            %    +(2*C.alpha-C.gamma+C.beta/2) CO+(C.omega/2) N2
            %                 <-V-> C_alpha H_beta O_gamma N_omega   [KPT_V]
            %
            % or
            %
            % C.alpha CO2+(C.gamma-2*C.alpha) H2O+(C.beta/2-C.gamma+2*C.alpha) H2
            % +(C.omega/2) N2 <-VI-> C_alpha H_beta O_gamma N_omega [KPT_VI]
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
            % low Hydrogen:    CO2 <-I -> CO+(1/2) O2              [KPT_I]
            % low Carbon:      H2O <-II-> H2+(1/2) O2             [KPT_II]
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Determination of the number of moles of the minor species from
            % the equilibrium condition for the above reaction
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if t
                if M.L_minor>0
                    Ni    = KPT_V .* NCO2.^(C.gamma-C.alpha-C.beta/2) .* NH2O.^(C.beta/2) ...
                        .* NCO.^(2*C.alpha-C.gamma+C.beta/2) .* NN2.^(C.omega/2) ...
                        .* DNfactor.^DNfactor_V;
                    Ni_old = N_IC_old(M.idx_minor,1)';
                    aux = find(Ni_old~=0);
                    if ~isempty(aux)
                        Ni(aux) = Ni_old(aux)+relax*(Ni(aux)-Ni_old(aux));
                    end
                    N_IC(M.idx_minor,1) = Ni;
                end
            
            NO2     = DNfactor*(KPT_I*NCO2/NCO)^2;
            NO2_old = N_IC_old(S.idx_O2,1);
            if NO2_old ~= 0
                NO2 = NO2_old+relax*(NO2-NO2_old);
            end
            
            % Correction of the number of moles of CO, H2, CO2, H2O and N2 from
            % atom conservation equations
            
            a = NCO2_0+NCO_0+NCgr_0-sum(N_IC(:,1).*A0(:,E.ind_C));
            b = NH2O_0+NH2_0-sum(N_IC(:,1).*A0(:,E.ind_H))/2;
            c = 2*NCO2_0+NCO_0+NH2O_0+2*NO2_0-2*NO2-sum(N_IC(:,1).*A0(:,E.ind_O));
            
            NCO2_old = NCO2;
            NCO_old  = NCO;
            NCgr_old = NCgr;
            NH2O_old = NH2O;
            NH2_old  = NH2;
            NN2_old  = NN2;
            
            mu = KPT_VII/DNfactor;
            
            NCgr = real(1/24*(24*a+(4*(2+KPT_IV))/(KPT_IV*mu)+(2*(1+1i*sqrt(3))*(-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu)))/(KPT_IV*(8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))+(2*1i*(1i+sqrt(3))*(8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))/(KPT_IV*mu^2)-(1/(6*KPT_IV^2*mu^3))*((2*(2+KPT_IV)*mu+((1+1i*sqrt(3))*mu^2*(-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu)))/(8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3)+1i*(1i+sqrt(3))*(8*mu^3-6*KPT_IV*mu^3-3*KPT_IV^2*mu^3+KPT_IV^3*mu^3-18*b*KPT_IV^2*mu^4-36*c*KPT_IV^2*mu^4-9*b*KPT_IV^3*mu^4+9*c*KPT_IV^3*mu^4+sqrt(mu^6*((-4+2*KPT_IV+KPT_IV^2*(-1+6*b*mu-6*c*mu))^3+(-8+6*KPT_IV+KPT_IV^3*(-1+9*b*mu-9*c*mu)+3*KPT_IV^2*(1+6*b*mu+12*c*mu))^2)))^(1/3))^2)));  
            NCO2 = (1+2*a*mu-2*mu*NCgr-sqrt(1+4*a*mu-4*mu*NCgr))/(2*mu);
            NCO = (-1+sqrt(1+4*a*mu-4*mu*NCgr))/(2*mu);
            NH2O = c-2*NCO2-NCO;
            NH2 = b-NH2O;
            
            NN2  = NN2_0-sum(N_IC(:,1).*A0(:,E.ind_N))/2; % N-atom conservation
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
            
            N_IC(S.idx_Cgr,1)= NCgr;
            N_IC(S.idx_CO,1) = NCO;
            N_IC(S.idx_CO2,1)= NCO2;
            N_IC(S.idx_H2O,1)= NH2O;
            N_IC(S.idx_H2,1) = NH2;
            N_IC(S.idx_O2,1) = NO2;
            N_IC(S.idx_N2,1) = NN2;
            N_IC(S.idx_He,1) = NHe;
            N_IC(S.idx_Ar,1) = NAr;
            
            % Overall number of moles in the product species
            
            NP = sum(N_IC(:,1).*(1-N_IC(:,10)));
%             NP = sum(N_IC(:,1));
            
            DeltaNP = norm([NP-NP_old, x-sum(N_IC(:,1).*A0(:,E.ind_C)), y-sum(N_IC(:,1).*A0(:,E.ind_H)), z-sum(N_IC(:,1).*A0(:,E.ind_O)), w-sum(N_IC(:,1).*A0(:,E.ind_N))]);
            
            if NCgr <= 1e-8
               N_IC = incomplete_phi_6(N_IC); 
               break; 
            end
        end
    end % TEST
end
