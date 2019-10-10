function [strThProp,E,S,M] = check_database(strMaster,strThProp,E,S,M)
M.Lminor = length(M.minor_products);
n_pass = [];
for n = 1:M.Lminor
    if ~any(strcmp(S.NameSpecies,M.minor_products{n}))
        n_pass = [n_pass, n];
    end
end
if ~isempty(n_pass)
    [E.Elements,E.NE] = set_elements(); % sets E.Elements list
    l_n_pass = length(n_pass);
    for i=1:l_n_pass
        Species = FullName2name(M.minor_products{n_pass(i)});
        
        if isfield(strMaster,Species)           
            % disp(SpeciesList{n_pass(i)})
            ctTInt = strMaster.(Species).ctTInt;
            tRange = strMaster.(Species).tRange;
            swtCondensed = sign(strMaster.(Species).swtCondensed);
            
            if ctTInt > 0
                
                [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = SpeciesThermProp(strMaster,M.minor_products{n_pass(i)},298.15,'molar',0);
                
                strThProp.(Species).name = Species;
                strThProp.(Species).FullName = M.minor_products{n_pass(i)};
                strThProp.(Species).txFormula = txFormula;
                strThProp.(Species).mm = mm;
                strThProp.(Species).hf = Hf0;
                strThProp.(Species).ef = Ef0;
                strThProp.(Species).swtCondensed = swtCondensed;
                
                T_vector   = [];
                DhT_vector = [];
                DeT_vector = [];
                s0_vector  = [];
                cp_vector  = [];
                cv_vector  = [];
                g0_vector  = [];
                
                Tmin = max(tRange{1}(1),200);
                Tmax = min(tRange{ctTInt}(2),6000);
                for T = [linspace(Tmin,298.15,4), linspace(350,Tmax,60)]
                    
                    [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = SpeciesThermProp(strMaster,M.minor_products{n_pass(i)},T,'molar',0);
                    T_vector   = [  T_vector; T     ];
                    DhT_vector = [DhT_vector; H0-Hf0];
                    DeT_vector = [DeT_vector; E0-Ef0];
                    s0_vector  = [ s0_vector; S0    ];
                    cp_vector  = [ cp_vector; Cp0   ];
                    cv_vector  = [ cv_vector; Cv0   ];
                    g0_vector  = [ g0_vector; DfG0  ];
                end
                
                strThProp.(Species).T   =   T_vector;
                strThProp.(Species).DhT = DhT_vector;
                strThProp.(Species).DeT = DeT_vector;
                strThProp.(Species).s0  =  s0_vector;
                strThProp.(Species).cp  =  cp_vector;
                strThProp.(Species).cv  =  cv_vector;
                strThProp.(Species).g0  =  g0_vector;
            else
                
                Tref = tRange(1);
                
                [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = SpeciesThermProp(strMaster,M.minor_products{n_pass(i)},Tref,'molar',0);
                
                strThProp.(Species).name = Species;
                strThProp.(Species).FullName = M.minor_products{n_pass(i)};
                strThProp.(Species).txFormula = txFormula;
                strThProp.(Species).mm  = mm;
                strThProp.(Species).hf  = Hf0;
                strThProp.(Species).ef  = Ef0;
                strThProp.(Species).T   = Tref;
                strThProp.(Species).DhT = [];
                strThProp.(Species).DeT = [];
                strThProp.(Species).s   = [];
                strThProp.(Species).cp  = [];
                strThProp.(Species).cv  = [];
                strThProp.(Species).g0  = [];
            end
            % INTERPOLATION CURVES
            strThProp.(Species).cPcurve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).cp);
            strThProp.(Species).cVcurve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).cv);
            strThProp.(Species).DeTcurve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).DeT);
            strThProp.(Species).DhTcurve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).DhT);
            strThProp.(Species).s0curve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).s0);
            strThProp.(Species).g0curve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).g0);
        else
            fprintf(['\n- Species ''',M.minor_products{n_pass(i)},''' does not exist as a field in strMaster structure ... '])
        end
    end
    
    % CONTAINED ELEMENTS
    S.NameSpecies = fieldnames(strThProp); % In case any of the minor products didnt exist in strMaster
    S.NSpecies = numel(S.NameSpecies);
    for k = length(S.NameSpecies):-1:1
        Species = S.NameSpecies{k};
        % Change uppercase 'L' to  lowercase 'l'
        Species(strfind(Species,'AL')+1)='l';
        Species(strfind(Species,'CL')+1)='l';
        Species(strfind(Species,'TL')+1)='l';
        Species(strfind(Species,'FL')+1)='l';
        % -----------------------------------------------
        Species(Species>='0' & Species<='9') = ' ';
        
        [idx0,idxf] = regexp(Species,"minus"); Species(idx0:idxf) = ' ';
        [idx0,idxf] = regexp(Species,"plus"); Species(idx0:idxf) = ' ';
        
        idx = find([(Species>='A' & Species<='Z'), true]);
        lgt = diff(idx);
        Tmp{k,1} = strtrim(mat2cell(Species, 1, lgt));
    end
    aux = unique(cat(2,Tmp{:}));
    n_pass = [];
    for n=length(aux):-1:1
        if any(strcmp(aux(n),E.Elements)) % Check Element existence
            n_pass = [n_pass, n];
        end
    end
    E.Elements = aux(n_pass);
    E.NE = numel(E.Elements);
    % INDEX OF EVALUABLE ELEMENTS
    E.ind_C = find(strcmp(E.Elements,'C'));
    E.ind_H = find(strcmp(E.Elements,'H'));
    E.ind_O = find(strcmp(E.Elements,'O'));
    E.ind_N = find(strcmp(E.Elements,'N'));
    E.ind_He = find(strcmp(E.Elements,'He'));
    E.ind_Ar = find(strcmp(E.Elements,'Ar'));
    % Element_matrix
    for i=S.NSpecies:-1:l_n_pass
        txFormula = strThProp.(S.NameSpecies{i,1}).txFormula;
        strThProp.(S.NameSpecies{i,1}).Element_matrix = set_element_matrix(txFormula);
        app.C.M0.Value(i,10) = strThProp.(app.S.NameSpecies{i,1}).swtCondensed;
    end
end

