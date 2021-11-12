function DB = generate_DB(DB_master)
    % Generate Database (DB) with thermochemical interpolation curves for
    % the species contained in DB_master
    if ~exist('DB', 'var')
        if exist('DB.mat', 'file')
            fprintf('NASA short database loaded from main path ... ')
            load('DB.mat', 'DB');
        else
            DB = get_DB(DB_master);
        end
        fprintf('OK!\n');
    else
        fprintf('NASA short database already loaded\n')
    end
end

% SUB-PASS FUNCTIONS
function DB = get_DB(DB_master)
    % Generate Database with thermochemical interpolation curves for the
    % species contained in DB_master
    LS = fieldnames(DB_master);

    fprintf('Generating short NASA database ... ')
    for i = 1:length(LS)
        species = FullName2name(LS{i});
        if isfield(DB_master,species)
            
            ctTInt = DB_master.(species).ctTInt;
            tRange = DB_master.(species).tRange;
            swtCondensed = sign(DB_master.(species).swtCondensed);
            
            if ctTInt > 0
                
                [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = get_speciesProperties(DB_master, LS{i}, 298.15, 'molar', 0);
                
                DB.(species).name = species;
                DB.(species).FullName = LS{i};
                DB.(species).txFormula = txFormula;
                DB.(species).mm = mm;
                DB.(species).hf = Hf0;
                DB.(species).ef = Ef0;
                DB.(species).swtCondensed = swtCondensed;
                
                NT   = 100;
                Tmin = max(tRange{1}(1), 200);
                Tmax = min(tRange{ctTInt}(2), 20000);
                T_vector = linspace(Tmin, Tmax, NT);
    
                for j = NT:-1:1
                    [~, ~, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, ~] = get_speciesProperties(DB_master, LS{i}, T_vector(j), 'molar', 0);
                    DhT_vector(j) = H0 - Hf0;
                    DeT_vector(j) = E0 - Ef0;
                    h0_vector(j)  = H0;
                    s0_vector(j)  = S0;
                    cp_vector(j)  = Cp0;
                    cv_vector(j)  = Cv0;
                    g0_vector(j)  = H0 - T_vector(j) * S0;
                end
                
                DB.(species).T = T_vector;
    
                % INTERPOLATION CURVES
                DB.(species).cPcurve  = griddedInterpolant(T_vector, cp_vector, 'pchip', 'linear');
                DB.(species).cVcurve  = griddedInterpolant(T_vector, cv_vector, 'pchip', 'linear');
                DB.(species).DhTcurve = griddedInterpolant(T_vector, DhT_vector, 'pchip', 'linear');
                DB.(species).DeTcurve = griddedInterpolant(T_vector, DeT_vector, 'pchip', 'linear');
                DB.(species).h0curve  = griddedInterpolant(T_vector, h0_vector, 'pchip', 'linear');
                DB.(species).s0curve  = griddedInterpolant(T_vector, s0_vector, 'pchip', 'linear');
                DB.(species).g0curve  = griddedInterpolant(T_vector, g0_vector, 'pchip', 'linear');
                
                % DATA COEFFICIENTS NASA 9 POLYNOMIAL (ONLY GASES)
                DB.(species).ctTInt = DB_master.(species).ctTInt;
                DB.(species).tRange = DB_master.(species).tRange;
                DB.(species).tExponents = DB_master.(species).tExponents;
                DB.(species).ctTInt = DB_master.(species).ctTInt;
                DB.(species).a = DB_master.(species).a;
                DB.(species).b  = DB_master.(species).b;
            else
                
                Tref = tRange(1);
                
                [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = get_speciesProperties(DB_master, LS{i}, Tref, 'molar', 0);
                
                DB.(species).name = species;
                DB.(species).FullName = LS{i};
                DB.(species).txFormula = txFormula;
                DB.(species).mm  = mm;
                DB.(species).hf  = Hf0;
                DB.(species).ef  = Ef0;
                DB.(species).swtCondensed = swtCondensed;
                DB.(species).T   = Tref;
                DB.(species).DhT = 0;
                DB.(species).DeT = 0;
                DB.(species).h0  = H0;
                DB.(species).s0  = S0;
                DB.(species).cp  = Cp0;
                DB.(species).cv  = Cv0;
                DB.(species).g0  = DfG0;
                DB.(species).ctTInt = 0;
            end
        else
            fprintf(['\n- Species ''', LS{i}, ''' does not exist as a field in strMaster structure ... '])
        end
    end
end