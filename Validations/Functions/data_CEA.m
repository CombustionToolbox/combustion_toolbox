function data = data_CEA(varargin)
    filename  = varargin{1};
    if ischar(filename)
        Nfiles = 1;
        filename = {filename};
    else
        Nfiles = length(filename);
    end
    data = [];
    j = 0;
    while j < Nfiles
%     for j = 1:length(Nfiles)
        clear mix
        j = j + 1;
        % READ DATA
        data_nasa = read_CEA(filename{j});
        % PROPERTIES
        mix.T = data_nasa.T;
        mix.rho = data_nasa.rho;
        mix.S = data_nasa.S;
        mix.g = data_nasa.G;
        mix.cP = data_nasa.cp; % cp [kJ/kg-K]
        mix.DhT = data_nasa.cp .* (data_nasa.T-298.15)*1e-3; % DhT_P [J/kg]
        if isfield(data_nasa,'rho2rho1')
            mix.u_preshock = data_nasa.u1;
            mix.u_postshock = data_nasa.u1 ./ data_nasa.rho2rho1;
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