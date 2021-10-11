function [Prop,phi,s] = data_CEA(varargin)
filename  = varargin{1};
data_nasa = read_CEA(filename);
% PROPERTIES
Prop.case(1,:) = data_nasa.T;
Prop.case(2,:) = data_nasa.rho;
Prop.case(3,:) = data_nasa.S;
Prop.case(4,:) = data_nasa.G;
Prop.case(5,:) = data_nasa.cp; % cp [kJ/kg-K]
Prop.case(6,:) = data_nasa.cp*(data_nasa.T-298.15)*1e-3; % DhT_P [J/kg]
Prop.case(6,:) = data_nasa.cp*(data_nasa.T-298.15)*1e-3; % DhT_P [J/kg]
if isfield(data_nasa,'rho2rho1')
    Prop.case(7,:) = data_nasa.u1;
    Prop.case(8,:) = data_nasa.u1./data_nasa.rho2rho1;
end
% EQUIVALENCE RATIO
phi = data_nasa.phi;
% MOLAR FRACTION SPECIES
if nargin > 1
    species = varargin{2};
    s.nasa = zeros(numel(data_nasa.X),numel(species));
    for i=numel(data_nasa.X):-1:1
        for t=numel(species):-1:1
            for k = length(data_nasa.X(i).mole):-1:1
                if strcmpi(data_nasa.X(i).mole{k,1},species{t})
                    s.nasa(i,t) = data_nasa.X(i).mole{k,2};
                end
            end
        end
    end
end
