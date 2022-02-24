function data = data_CEA(varargin)
    filename  = varargin{1};
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
    while j < Nfiles
%     for j = 1:length(Nfiles)
        clear mix
        j = j + 1;
        % READ DATA
        data_nasa = read_CEA(filename{j});
        try
            % PROPERTIES
            mix.p = data_nasa.P; % [bar]
            mix.T = data_nasa.T; % [K]
            mix.rho = data_nasa.rho; % [kg/m3]
            mix.h = data_nasa.H; % [kJ/kg]
            mix.e = data_nasa.U; % [kJ/kg]
            mix.S = data_nasa.S; % [kJ/kg-K]
            mix.g = data_nasa.G; % [kJ/kg]
            mix.cP = data_nasa.cp; % [kJ/kg-K]
            mix.gamma_s = data_nasa.gamma_s; % [-]
            mix.cV = data_nasa.cp ./ mix.gamma_s; % [kJ/kg-K]
            mix.DhT = data_nasa.cp .* (data_nasa.T-298.15); % [kJ/kg]
            if isfield(data_nasa,'rho2rho1')
                mix.u_preshock = data_nasa.u1; % [m/s]
                mix.u_postshock = data_nasa.u1 ./ data_nasa.rho2rho1; % [m/s]
            end
            % EQUIVALENCE RATIO
            mix.phi = data_nasa.phi;
            % MOLAR FRACTION SPECIES
            if nargin > 1
                species = varargin{2};
                NS = length(species);
                mix.Xi = zeros(NS, length(data_nasa.X));
                for i = 1:length(data_nasa.X)
                    for t = 1:NS
                        for k = length(data_nasa.X(i).mole):-1:1
                            if strcmpi(data_nasa.X(i).mole{k, 1}, species{t})
                                mix.Xi(t, i) = data_nasa.X(i).mole{k, 2};
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
            mix1.cV = data_nasa.cp1 ./ mix1.gamma_s; % [kJ/kg-K]
            mix1.DhT = data_nasa.cp1 .* (data_nasa.T1 - 298.15); % [kJ/kg]
            mix1.u = data_nasa.u1; % [m/s]
            % EQUIVALENCE RATIO
            mix1.phi = data_nasa.phi;

            % PROPERTIES MIX 2
            mix2.p = data_nasa.P2; % [bar]
            mix2.T = data_nasa.T2; % [K]
            mix2.rho = data_nasa.rho2; % [kg/m3]
            mix2.h = data_nasa.H2; % [kJ/kg]
            mix2.e = data_nasa.U2; % [kJ/kg]
            mix2.S = data_nasa.S2; % [kJ/kg-K]
            mix2.g = data_nasa.G2; % [kJ/kg]
            mix2.cP = data_nasa.cp2; % [kJ/kg-K]
            mix2.gamma_s = data_nasa.gamma_s2; % [-]
            mix2.cV = data_nasa.cp2 ./ mix2.gamma_s; % [kJ/kg-K]
            mix2.DhT = data_nasa.cp2 .* (data_nasa.T2 - 298.15); % [kJ/kg]
            mix2.u = data_nasa.u2; % [m/s]
            mix2.dVdp_T = data_nasa.dVdp_T; % [-]
            mix2.dVdT_p = data_nasa.dVdT_p; % [-]
            % EQUIVALENCE RATIO
            mix2.phi = data_nasa.phi;
            % MOLAR FRACTION SPECIES MIX 2
            if nargin > 1
                species = varargin{2};
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
    if isempty(data1)
        data1 = data2;
        return
    end
    varname = fieldnames(data1);
    N = length(varname);
    for i = 1:N
        data1.(varname{i}) = [data1.(varname{i}), data2.(varname{i})];
    end
end