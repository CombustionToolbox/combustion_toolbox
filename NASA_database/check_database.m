function [strThProp,E,S,C] = check_database(varargin)
%% CHECK THAT SPECIES ARE IN DB 

if nargin == 3
    self = varargin{1}; strMaster = varargin{2}; strThProp = varargin{3};
    E = self.E; S = self.S; C = self.C;

    n_pass = [];
    for n = 1:S.NS
        if ~any(strcmp(S.LS_DB, S.LS{n}))
            n_pass = [n_pass, n];
        end
    end
    if ~isempty(n_pass)
        [E.elements,E.NE] = set_elements(); % sets E.elements list
        l_n_pass = length(n_pass);
        for i=1:l_n_pass
            Species = FullName2name(S.LS{n_pass(i)});
            
            if isfield(strMaster, Species)
                % disp(SpeciesList{n_pass(i)})
                ctTInt = strMaster.(Species).ctTInt;
                tRange = strMaster.(Species).tRange;
                swtCondensed = sign(strMaster.(Species).swtCondensed);
                
                if ctTInt > 0
                    
                    [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = SpeciesThermProp(strMaster, S.LS{n_pass(i)},298.15,'molar',0);
            
                    strThProp.(Species).name = Species;
                    strThProp.(Species).FullName = S.LS{n_pass(i)};
                    strThProp.(Species).txFormula = txFormula;
                    strThProp.(Species).mm = mm;
                    strThProp.(Species).hf = Hf0;
                    strThProp.(Species).ef = Ef0;
                    strThProp.(Species).swtCondensed = swtCondensed;

                    T_vector   = [];
                    DhT_vector = [];
                    DeT_vector = [];
                    h0_vector  = [];
                    s0_vector  = [];
                    cp_vector  = [];
                    cv_vector  = [];
                    g0_vector  = [];

                    Tmin = max(tRange{1}(1), 200);
                    Tmax = min(tRange{ctTInt}(2), 20000);
                    if abs(Tmin - 298.15) < 1e4 
                        Trange1 = 298.15;
                    else
                        Trange1 = linspace(Tmin, 298.15, 10);
                    end
                    Trange2 = linspace(350, Tmax, 100);
                    for T = [Trange1, Trange2]
                        [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = SpeciesThermProp(strMaster, S.LS{i},T,'molar',0);
                        T_vector   = [  T_vector; T     ];
                        DhT_vector = [DhT_vector; H0-Hf0];
                        DeT_vector = [DeT_vector; E0-Ef0];
                        h0_vector  = [h0_vector;  H0    ];
                        s0_vector  = [s0_vector;  S0    ];
                        cp_vector  = [cp_vector;  Cp0   ];
                        cv_vector  = [cv_vector;  Cv0   ];
                        g0_vector  = [g0_vector;  H0 - T*S0];
                    end

                    strThProp.(Species).T   = T_vector;
                    strThProp.(Species).DhT = DhT_vector; 
                    strThProp.(Species).DeT = DeT_vector;
                    strThProp.(Species).h0  = h0_vector; 
                    strThProp.(Species).s0  = s0_vector;
                    strThProp.(Species).cp  = cp_vector;
                    strThProp.(Species).cv  = cv_vector;
                    strThProp.(Species).g0  = g0_vector;

                    % INTERPOLATION CURVES
                    strThProp.(Species).cPcurve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).cp, 'pchip', 'pchip');
                    strThProp.(Species).cVcurve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).cv, 'pchip', 'pchip');
                    strThProp.(Species).DhTcurve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).DhT, 'pchip', 'pchip');
                    strThProp.(Species).DeTcurve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).DeT, 'pchip', 'pchip');
                    strThProp.(Species).h0curve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).h0, 'pchip', 'pchip');
                    strThProp.(Species).s0curve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).s0, 'pchip', 'pchip');
                    strThProp.(Species).g0curve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).g0, 'pchip', 'pchip');

                    % DATA COEFFICIENTS NASA 9 POLYNOMIAL
                    strThProp.(Species).ctTInt = strMaster.(Species).ctTInt;
                    strThProp.(Species).tRange = strMaster.(Species).tRange;
                    strThProp.(Species).tExponents = strMaster.(Species).tExponents;
                    strThProp.(Species).ctTInt = strMaster.(Species).ctTInt;
                    strThProp.(Species).a = strMaster.(Species).a;
                    strThProp.(Species).b  = strMaster.(Species).b;
                else
                    
                    Tref = tRange(1);
                    
                    [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = SpeciesThermProp(strMaster,S.LS{n_pass(i)},Tref,'molar',0);
                    
                    strThProp.(Species).name = Species;
                    strThProp.(Species).FullName = S.LS{n_pass(i)};
                    strThProp.(Species).txFormula = txFormula;
                    strThProp.(Species).mm  = mm;
                    strThProp.(Species).hf  = Hf0;
                    strThProp.(Species).ef  = Ef0;
                    strThProp.(Species).swtCondensed = swtCondensed;
                    strThProp.(Species).T   = Tref;
                    strThProp.(Species).DhT = [];
                    strThProp.(Species).DeT = [];
                    strThProp.(Species).h0  = [];
                    strThProp.(Species).s   = [];
                    strThProp.(Species).cp  = [];
                    strThProp.(Species).cv  = [];
                    strThProp.(Species).g0  = [];
                end
            else
                fprintf(['\n- Species ''', S.LS{n_pass(i)},''' does not exist as a field in strMaster structure ... '])
            end
        end
        
        % CONTAINED ELEMENTS
        S.LS_DB = fieldnames(strThProp); % In case any of the minor products didnt exist in strMaster
        S.NS_DB = numel(S.LS_DB);
        for k = S.NS_DB:-1:1
            Species = S.LS_DB{k};
            % Change uppercase 'L' to  lowercase 'l'
            Species(strfind(Species,'AL')+1)='l';
            Species(strfind(Species,'CL')+1)='l';
            Species(strfind(Species,'TL')+1)='l';
            Species(strfind(Species,'FL')+1)='l';
            % -----------------------------------------------
            Species(strfind(Species,'eminus')) = 'E';
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
            if any(strcmp(aux(n),E.elements)) % Check Element existence
                n_pass = [n_pass, n];
            end
        end
        E.elements = aux(n_pass);
        E.NE = numel(E.elements);
        % INDEX OF EVALUABLE ELEMENTS
        E.ind_C = find(strcmp(E.elements,'C'));
        E.ind_H = find(strcmp(E.elements,'H'));
        E.ind_O = find(strcmp(E.elements,'O'));
        E.ind_N = find(strcmp(E.elements,'N'));
        E.ind_He = find(strcmp(E.elements,'He'));
        E.ind_Ar = find(strcmp(E.elements,'Ar'));
        E.ind_E = find(strcmp(E.elements,'E'));
        % Element_matrix
        for i=S.NS:-1:l_n_pass
            txFormula = strThProp.(S.LS_DB{i,1}).txFormula;
            strThProp.(S.LS_DB{i,1}).Element_matrix = set_element_matrix(txFormula,E.elements);
            C.M0.Value(i,10) = strThProp.(S.LS_DB{i,1}).swtCondensed;
        end
    end
    %%%% CHECK SPECIFIC LIST OF SPECIES
elseif nargin == 4
    self = varargin{1}; strMaster = varargin{2}; strThProp = varargin{3};
    E = self.E; S = self.S; C = self.C; LS_check = varargin{4};
    
    Lminors = length(LS_check);
    n_pass = [];
    for n = 1:Lminors
        if ~any(strcmp(S.LS_DB, LS_check{n}))
            n_pass = [n_pass, n];
        end
    end
    if ~isempty(n_pass)
        [E.elements,E.NE] = set_elements(); % sets E.elements list
        l_n_pass = length(n_pass);
        for i=1:l_n_pass
            Species = FullName2name(LS_check{n_pass(i)});
            
            if isfield(strMaster,Species)
                % disp(SpeciesList{n_pass(i)})
                ctTInt = strMaster.(Species).ctTInt;
                tRange = strMaster.(Species).tRange;
                swtCondensed = sign(strMaster.(Species).swtCondensed);
                
                if ctTInt > 0
                    
                    [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = SpeciesThermProp(strMaster,LS_check{n_pass(i)},298.15,'molar',0);
                    
                    strThProp.(Species).name = Species;
                    strThProp.(Species).FullName = LS_check{n_pass(i)};
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
                    Tmax = min(tRange{ctTInt}(2), 20000);
                    for T = [linspace(Tmin,298.15,4), linspace(350,Tmax,60)]
                        
                        [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = SpeciesThermProp(strMaster,LS_check{n_pass(i)},T,'molar',0);
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
                    
                    [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = SpeciesThermProp(strMaster,LS_check{n_pass(i)},Tref,'molar',0);
                    
                    strThProp.(Species).name = Species;
                    strThProp.(Species).FullName = LS_check{n_pass(i)};
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
                fprintf(['\n- Species ''',S.LS{n_pass(i)},''' does not exist as a field in strMaster structure ... '])
            end
        end
        
        % CONTAINED ELEMENTS
        S.LS_DB = fieldnames(strThProp); % In case any of the minor products didnt exist in strMaster
        S.NS = numel(S.LS_DB);
        for k = length(S.LS_DB):-1:1
            Species = S.LS_DB{k};
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
            if any(strcmp(aux(n),E.elements)) % Check Element existence
                n_pass = [n_pass, n];
            end
        end
        E.elements = aux(n_pass);
        E.NE = numel(E.elements);
        % INDEX OF EVALUABLE ELEMENTS
        E.ind_C = find(strcmp(E.elements,'C'));
        E.ind_H = find(strcmp(E.elements,'H'));
        E.ind_O = find(strcmp(E.elements,'O'));
        E.ind_N = find(strcmp(E.elements,'N'));
        E.ind_He = find(strcmp(E.elements,'He'));
        E.ind_Ar = find(strcmp(E.elements,'Ar'));
        % Element_matrix
        for i=S.NS:-1:l_n_pass
            txFormula = strThProp.(S.LS_DB{i,1}).txFormula;
            strThProp.(S.LS_DB{i,1}).Element_matrix = set_element_matrix(txFormula,E.elements);
            C.M0.Value(i,10) = strThProp.(S.LS_DB{i,1}).swtCondensed;
        end
    end
end