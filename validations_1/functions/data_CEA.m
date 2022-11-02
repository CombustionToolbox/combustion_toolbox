function data = data_CEA(filename, varargin)

    if ischar(filename)
        Nfiles = 1;
        filename = {filename};
    else
        Nfiles = length(filename);
    end

    data = [];
    data1 = []; % Only for shocks
    data2 = []; % Only for shocks
    j = 0;
    l = 0;

    while j < Nfiles
        %     for j = 1:length(Nfiles)
        clear mix
        j = j + 1;
        % READ DATA
        data_nasa = read_CEA(filename{j});

        try
            % data_nasa = append_fields_struct(data_nasa);
            % PROPERTIES
            mix.p = data_nasa.P(:, end - l)'; % [bar]
            mix.T = data_nasa.T(:, end - l)'; % [K]
            mix.rho = data_nasa.rho(:, end - l)'; % [kg/m3]
            mix.h = data_nasa.H(:, end - l)'; % [kJ/kg]
            mix.e = data_nasa.U(:, end - l)'; % [kJ/kg]
            mix.S = data_nasa.S(:, end - l)'; % [kJ/kg-K]
            mix.g = data_nasa.G(:, end - l)'; % [kJ/kg]
            mix.cP = data_nasa.cp(:, end - l)'; % [kJ/kg-K]
            mix.gamma_s = data_nasa.gamma_s(:, end - l)'; % [-]
            mix.cV = data_nasa.cp(:, end - l)' ./ mix.gamma_s; % [kJ/kg-K]
            mix.DhT = data_nasa.cp(:, end - l)' .* (data_nasa.T(:, end - l)' - 298.15); % [kJ/kg]
            mix.W = data_nasa.W(:, end - l)'; % [g/mol]
            mix.sound = data_nasa.sound(:, end - l)'; % [m/s]

            if isfield(data_nasa, 'dVdp_T')
                mix.dVdp_T = data_nasa.dVdp_T(:, end - l)'; % [-]
                mix.dVdT_p = data_nasa.dVdT_p(:, end - l)'; % [-]
            end

            if isfield(data_nasa, 'rho2rho1')

                try % reflected
                    mix.u_postshock = data_nasa.u2 - data_nasa.v2; % [m/s]
                    mix.u_preshock = data_nasa.u1; % [m/s]
                    mix.u = mix.u_preshock; % [m/s]
                catch % incident
                    mix.u_preshock = data_nasa.u1; % [m/s]
                    mix.u_postshock = data_nasa.u1 ./ data_nasa.rho2rho1; % [m/s]
                    mix.u = mix.u_preshock; % [m/s]
                end

            end

            if isfield(data_nasa, 'Aratio')
                mix.Aratio = data_nasa.Aratio(:, end - l)'; % [-]
                mix.cstar = data_nasa.cstar(:, end - l)'; % [m/s]
                mix.cf = data_nasa.cf(:, end - l)'; % [-]
                mix.I_vac = data_nasa.I_vac(:, end - l)' / 9.80665; % [s]
                mix.I_sp = data_nasa.I_sp(:, end - l)' / 9.80665; % [s]
                mix.u = data_nasa.u(:, end - l)'; % [m/s]
            end

            % EQUIVALENCE RATIO
            mix.phi = data_nasa.phi; % [-]
            % MOLAR FRACTION SPECIES
            if nargin > 1
                species = varargin{1};
                NS = length(species);

                if ~isstruct(data_nasa.X)
                    index = find_ind(data_nasa.LS, species);
                    mix.Xi = data_nasa.X(index, :);
                else
                    mix.Xi = zeros(NS, length(data_nasa.X));

                    for i = 1:length(data_nasa.X)

                        for t = 1:NS

                            for k = length(data_nasa.X(i).mole):-1:1

                                if strcmpi(data_nasa.X(i).mole{k, 1}, species{t})
                                    mix.Xi(t, i) = data_nasa.X(i).mole{k, 2}(end);
                                end

                            end

                        end

                    end

                end

            end

            data = join_datacell(data, mix);
        catch
            % PROPERTIES MIX 1
            mix1.p = data_nasa.P1; % [bar]
            mix1.T = data_nasa.T1; % [K]
            mix1.rho = data_nasa.rho1; % [kg/m3]
            mix1.h = data_nasa.H1; % [kJ/kg]
            mix1.e = data_nasa.U1; % [kJ/kg]
            mix1.S = data_nasa.S1; % [kJ/kg-K]
            mix1.g = data_nasa.G1; % [kJ/kg]
            mix1.cP = data_nasa.cp1; % [kJ/kg-K]
            mix1.gamma_s = data_nasa.gamma_s1; % [-]
            mix1.sound = data_nasa.sound1; % [m/s]
            mix1.cV = data_nasa.cp1 ./ mix1.gamma_s; % [kJ/kg-K]
            mix1.DhT = data_nasa.cp1 .* (data_nasa.T1 - 298.15); % [kJ/kg]
            mix1.u = data_nasa.u1; % [m/s]
            mix1.u_preshock = mix1.u; % [m/s]
            mix1.W = data_nasa.W1; % [g/mol]
            % EQUIVALENCE RATIO
            mix1.phi = data_nasa.phi; % [-]

            % PROPERTIES MIX 2
            mix2.p = data_nasa.P2; % [bar]
            mix2.T = data_nasa.T2; % [K]
            mix2.rho = data_nasa.rho2; % [kg/m3]

            if length(mix2.rho) ~= length(mix2.T)
                mix2.rho = mix2.rho(mix2.rho ~= 1);
            end

            mix2.h = data_nasa.H2; % [kJ/kg]
            mix2.e = data_nasa.U2; % [kJ/kg]
            mix2.S = data_nasa.S2; % [kJ/kg-K]
            mix2.g = data_nasa.G2; % [kJ/kg]
            mix2.cP = data_nasa.cp2; % [kJ/kg-K]
            mix2.gamma_s = data_nasa.gamma_s2; % [-]
            mix2.sound = data_nasa.sound2; % [m/s]
            mix2.cV = data_nasa.cp2 ./ mix2.gamma_s; % [kJ/kg-K]
            mix2.DhT = data_nasa.cp2 .* (data_nasa.T2 - 298.15); % [kJ/kg]
            mix2.u = data_nasa.u2; % [m/s]
            mix2.u_preshock = mix1.u; % [m/s]
            mix2.W = data_nasa.W2; % [g/mol]

            if isfield(data_nasa, 'v_shock')
                mix2.v_shock = data_nasa.v_shock; % [m/s]
                mix2.u_postshock = mix2.v_shock; % [m/s]
                mix1.u_postshock = mix2.v_shock; % [m/s]
            end

            mix2.dVdp_T = data_nasa.dVdp_T; % [-]
            mix2.dVdT_p = data_nasa.dVdT_p; % [-]
            % EQUIVALENCE RATIO
            mix2.phi = data_nasa.phi; % [-]
            % MOLAR FRACTION SPECIES MIX 2
            if nargin > 1
                species = varargin{1};
                NS = length(species);
                mix2.Xi = zeros(NS, length(data_nasa.X));

                for t = 1:NS

                    for k = length(data_nasa.X(1).mole):-1:1

                        for i = length(data_nasa.X(1).mole{k, 2}):-1:1

                            if strcmpi(data_nasa.X(1).mole{k, 1}, species{t})
                                mix2.Xi(t, i) = data_nasa.X(1).mole{k, 2}(i);
                            end

                        end

                    end

                end

            end

            data1 = join_datacell(data1, mix1);
            data2 = join_datacell(data2, mix2);
        end

    end

    if ~isempty(data1)
        data.mix1 = data1;
        data.mix2 = data2;
    end

end

% SUB-PASS FUNCTIONS
function data1 = join_datacell(data1, data2)

    try

        if isempty(data1)
            data1 = data2;
            return
        end

        varname = fieldnames(data1);
        N = length(varname);

        for i = 1:N
            data1.(varname{i}) = [data1.(varname{i}), data2.(varname{i})];
        end

    catch
        disp('error')
    end

end
