function [DB, E, S, C] = check_DB(varargin)
% CHECK THAT SPECIES ARE IN DATABASE 

if nargin == 3
    self = varargin{1}; DB_master = varargin{2}; DB = varargin{3};
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
            
            if isfield(DB_master, Species)
                % disp(SpeciesList{n_pass(i)})
                ctTInt = DB_master.(Species).ctTInt;
                tRange = DB_master.(Species).tRange;
                swtCondensed = sign(DB_master.(Species).swtCondensed);
                
                if ctTInt > 0
                    
                    [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = get_speciesProperties(DB_master, S.LS{n_pass(i)},298.15,'molar',0);
            
                    DB.(Species).name = Species;
                    DB.(Species).FullName = S.LS{n_pass(i)};
                    DB.(Species).txFormula = txFormula;
                    DB.(Species).mm = mm;
                    DB.(Species).hf = Hf0;
                    DB.(Species).ef = Ef0;
                    DB.(Species).swtCondensed = swtCondensed;

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
                        [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = get_speciesProperties(DB_master, S.LS{i},T,'molar',0);
                        T_vector   = [  T_vector; T     ];
                        DhT_vector = [DhT_vector; H0-Hf0];
                        DeT_vector = [DeT_vector; E0-Ef0];
                        h0_vector  = [h0_vector;  H0    ];
                        s0_vector  = [s0_vector;  S0    ];
                        cp_vector  = [cp_vector;  Cp0   ];
                        cv_vector  = [cv_vector;  Cv0   ];
                        g0_vector  = [g0_vector;  H0 - T*S0];
                    end

                    DB.(Species).T   = T_vector;
                    DB.(Species).DhT = DhT_vector; 
                    DB.(Species).DeT = DeT_vector;
                    DB.(Species).h0  = h0_vector; 
                    DB.(Species).s0  = s0_vector;
                    DB.(Species).cp  = cp_vector;
                    DB.(Species).cv  = cv_vector;
                    DB.(Species).g0  = g0_vector;

                    % INTERPOLATION CURVES
                    DB.(Species).cPcurve = griddedInterpolant(DB.(Species).T,DB.(Species).cp, 'pchip', 'linear');
                    DB.(Species).cVcurve = griddedInterpolant(DB.(Species).T,DB.(Species).cv, 'pchip', 'linear');
                    DB.(Species).DhTcurve = griddedInterpolant(DB.(Species).T,DB.(Species).DhT, 'pchip', 'linear');
                    DB.(Species).DeTcurve = griddedInterpolant(DB.(Species).T,DB.(Species).DeT, 'pchip', 'linear');
                    DB.(Species).h0curve = griddedInterpolant(DB.(Species).T,DB.(Species).h0, 'pchip', 'linear');
                    DB.(Species).s0curve = griddedInterpolant(DB.(Species).T,DB.(Species).s0, 'pchip', 'linear');
                    DB.(Species).g0curve = griddedInterpolant(DB.(Species).T,DB.(Species).g0, 'pchip', 'linear');

                    % DATA COEFFICIENTS NASA 9 POLYNOMIAL
                    DB.(Species).ctTInt = DB_master.(Species).ctTInt;
                    DB.(Species).tRange = DB_master.(Species).tRange;
                    DB.(Species).tExponents = DB_master.(Species).tExponents;
                    DB.(Species).ctTInt = DB_master.(Species).ctTInt;
                    DB.(Species).a = DB_master.(Species).a;
                    DB.(Species).b  = DB_master.(Species).b;
                else
                    
                    Tref = tRange(1);
                    
                    [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = get_speciesProperties(DB_master,S.LS{n_pass(i)},Tref,'molar',0);
                    
                    DB.(Species).name = Species;
                    DB.(Species).FullName = S.LS{n_pass(i)};
                    DB.(Species).txFormula = txFormula;
                    DB.(Species).mm  = mm;
                    DB.(Species).hf  = Hf0;
                    DB.(Species).ef  = Ef0;
                    DB.(Species).swtCondensed = swtCondensed;
                    DB.(Species).T   = Tref;
                    DB.(Species).DhT = [];
                    DB.(Species).DeT = [];
                    DB.(Species).h0  = [];
                    DB.(Species).s   = [];
                    DB.(Species).cp  = [];
                    DB.(Species).cv  = [];
                    DB.(Species).g0  = [];
                end
            else
                fprintf(['\n- Species ''', S.LS{n_pass(i)},''' does not exist as a field in DB_master structure ... '])
            end
        end
        
        % CONTAINED ELEMENTS
        S.LS_DB = fieldnames(DB); % In case any of the minor products didnt exist in DB_master
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
            txFormula = DB.(S.LS_DB{i,1}).txFormula;
            DB.(S.LS_DB{i,1}).Element_matrix = set_element_matrix(txFormula,E.elements);
            C.M0.Value(i,10) = DB.(S.LS_DB{i,1}).swtCondensed;
        end
    end
    %%%% CHECK SPECIFIC LIST OF SPECIES
elseif nargin == 4
    self = varargin{1}; DB_master = varargin{2}; DB = varargin{3};
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
            
            if isfield(DB_master,Species)
                % disp(SpeciesList{n_pass(i)})
                ctTInt = DB_master.(Species).ctTInt;
                tRange = DB_master.(Species).tRange;
                swtCondensed = sign(DB_master.(Species).swtCondensed);
                
                if ctTInt > 0
                    
                    [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = get_speciesProperties(DB_master,LS_check{n_pass(i)},298.15,'molar',0);
                    
                    DB.(Species).name = Species;
                    DB.(Species).FullName = LS_check{n_pass(i)};
                    DB.(Species).txFormula = txFormula;
                    DB.(Species).mm = mm;
                    DB.(Species).hf = Hf0;
                    DB.(Species).ef = Ef0;
                    DB.(Species).swtCondensed = swtCondensed;
                    
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
                        
                        [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = get_speciesProperties(DB_master,LS_check{n_pass(i)},T,'molar',0);
                        T_vector   = [  T_vector; T     ];
                        DhT_vector = [DhT_vector; H0-Hf0];
                        DeT_vector = [DeT_vector; E0-Ef0];
                        s0_vector  = [ s0_vector; S0    ];
                        cp_vector  = [ cp_vector; Cp0   ];
                        cv_vector  = [ cv_vector; Cv0   ];
                        g0_vector  = [ g0_vector; DfG0  ];
                    end
                    
                    DB.(Species).T   =   T_vector;
                    DB.(Species).DhT = DhT_vector;
                    DB.(Species).DeT = DeT_vector;
                    DB.(Species).s0  =  s0_vector;
                    DB.(Species).cp  =  cp_vector;
                    DB.(Species).cv  =  cv_vector;
                    DB.(Species).g0  =  g0_vector;
                else
                    
                    Tref = tRange(1);
                    
                    [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = get_speciesProperties(DB_master,LS_check{n_pass(i)},Tref,'molar',0);
                    
                    DB.(Species).name = Species;
                    DB.(Species).FullName = LS_check{n_pass(i)};
                    DB.(Species).txFormula = txFormula;
                    DB.(Species).mm  = mm;
                    DB.(Species).hf  = Hf0;
                    DB.(Species).ef  = Ef0;
                    DB.(Species).T   = Tref;
                    DB.(Species).DhT = [];
                    DB.(Species).DeT = [];
                    DB.(Species).s   = [];
                    DB.(Species).cp  = [];
                    DB.(Species).cv  = [];
                    DB.(Species).g0  = [];
                end
                % INTERPOLATION CURVES
                DB.(Species).cPcurve = griddedInterpolant(DB.(Species).T,DB.(Species).cp);
                DB.(Species).cVcurve = griddedInterpolant(DB.(Species).T,DB.(Species).cv);
                DB.(Species).DeTcurve = griddedInterpolant(DB.(Species).T,DB.(Species).DeT);
                DB.(Species).DhTcurve = griddedInterpolant(DB.(Species).T,DB.(Species).DhT);
                DB.(Species).s0curve = griddedInterpolant(DB.(Species).T,DB.(Species).s0);
                DB.(Species).g0curve = griddedInterpolant(DB.(Species).T,DB.(Species).g0);
            else
                fprintf(['\n- Species ''',S.LS{n_pass(i)},''' does not exist as a field in DB_master structure ... '])
            end
        end
        
        % CONTAINED ELEMENTS
        S.LS_DB = fieldnames(DB); % In case any of the minor products didnt exist in DB_master
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
            txFormula = DB.(S.LS_DB{i,1}).txFormula;
            DB.(S.LS_DB{i,1}).Element_matrix = set_element_matrix(txFormula,E.elements);
            C.M0.Value(i,10) = DB.(S.LS_DB{i,1}).swtCondensed;
        end
    end
end