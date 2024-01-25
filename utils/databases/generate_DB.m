function DB = generate_DB(DB_master, varargin)
    % Generate Database (DB) with thermochemical interpolation curves for
    % the species contained in DB_master
    %
    % Args:
    %     DB_master (struct): Database with the thermodynamic data of the chemical species
    %
    % Optional Args:
    %     LS (cell): List of species to be included in DB
    %
    % Returns:
    %     DB (struct): Database with custom thermodynamic polynomials functions generated from NASAs 9 polynomials fits
    %
    % Examples:
    %     * DB = generate_DB(DB_master)
    %     * DB = generate_DB(DB_master, {'CO2', 'H2O', 'O2', 'N2'})
    
    % Default
    LS = fieldnames(DB_master);

    % Check initial inputs
    if nargin > 1
        assert(iscell(varargin{1}), 'List of species must be a cell.');
            
        LS = varargin{1};
    end

    % Load database
    if exist('DB.mat', 'file') && nargin == 1
        fprintf('NASA database with thermo loaded from main path ... ')
        load('DB.mat', 'DB');
    else
        DB = get_DB(DB_master, LS);
    end

    fprintf('OK!\n');
end

% SUB-PASS FUNCTIONS
function DB = get_DB(DB_master, LS)
    % Generate Database (DB) with thermochemical interpolation curves for
    % the species contained in DB_master

    fprintf('Generating NASA database with thermo ... ')

    for i = 1:length(LS)
        species = FullName2name(LS{i});

        if isfield(DB_master, species)
            ctTInt = DB_master.(species).ctTInt;
            tRange = DB_master.(species).tRange;
            phase = sign(DB_master.(species).phase);
            
            % Species that can be evaluated at different temperatures
            if ctTInt > 0

                [txFormula, mm, ~, Hf0, ~, Ef0, ~, ~] = get_speciesProperties(DB_master, LS{i}, 298.15, 'molar', 0);

                DB.(species).FullName = DB_master.(species).FullName;
                DB.(species).name = species;
                DB.(species).comments = DB_master.(species).comments;
                DB.(species).txFormula = txFormula;
                DB.(species).mm = mm;
                DB.(species).hf = Hf0;
                DB.(species).ef = Ef0;
                DB.(species).phase = phase;

                NT = 200;
                Tmin = tRange{1}(1);
                Tmax = tRange{ctTInt}(2);
                T_vector = linspace(Tmin, Tmax, NT);

                for j = NT:-1:1
                    [~, ~, Cp0, ~, H0, ~, S0, ~] = get_speciesProperties(DB_master, LS{i}, T_vector(j), 'molar', 0);
                    cp_vector(j) = Cp0;
                    h0_vector(j) = H0;
                    s0_vector(j) = S0;
                    g0_vector(j) = H0 - T_vector(j) * S0;
                end

                DB.(species).T = T_vector;

                % Interpolation curves
                DB.(species).cPcurve = griddedInterpolant(T_vector, cp_vector, 'pchip', 'linear');
                DB.(species).h0curve = griddedInterpolant(T_vector, h0_vector, 'pchip', 'linear');
                DB.(species).s0curve = griddedInterpolant(T_vector, s0_vector, 'pchip', 'linear');
                DB.(species).g0curve = griddedInterpolant(T_vector, g0_vector, 'pchip', 'linear');

                % Coefficients NASA's 9 polynomial fits
                DB.(species).ctTInt = DB_master.(species).ctTInt;
                DB.(species).tRange = DB_master.(species).tRange;
                DB.(species).tExponents = DB_master.(species).tExponents;
                DB.(species).ctTInt = DB_master.(species).ctTInt;
                DB.(species).a = DB_master.(species).a;
                DB.(species).b = DB_master.(species).b;

            % Species that can be evaluated at a fixed temperature
            else

                Tref = tRange(1);

                [txFormula, mm, Cp0, Hf0, H0, Ef0, S0, DfG0] = get_speciesProperties(DB_master, LS{i}, Tref, 'molar', 0);

                DB.(species).FullName = DB_master.(species).FullName;
                DB.(species).name = species;
                DB.(species).comments = DB_master.(species).comments;
                DB.(species).txFormula = txFormula;
                DB.(species).mm = mm;
                DB.(species).hf = Hf0;
                DB.(species).ef = Ef0;
                DB.(species).phase = phase;
                DB.(species).T = Tref;
                DB.(species).h0 = H0;
                DB.(species).s0 = S0;
                DB.(species).cp = Cp0;
                DB.(species).g0 = DfG0;
                DB.(species).ctTInt = 0;
            end

        else
            fprintf(['\n- Species ''', LS{i}, ''' does not exist as a field in DB_master structure ... '])
        end

    end

end
