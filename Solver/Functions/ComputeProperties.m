function mix =  ComputeProperties(self, SpeciesMatrix, p, T)

R0 = self.C.R0; % [J/(K mol)]. Universal gas constant
mix.NatomE = sum(SpeciesMatrix(:, 1) .* self.C.A0.value);

if isempty(self.E.ind_C), mix.x = 0; else, mix.x = mix.NatomE(self.E.ind_C); end
if isempty(self.E.ind_H), mix.y = 0; else, mix.y = mix.NatomE(self.E.ind_H); end
if isempty(self.E.ind_O), mix.z = 0; else, mix.z = mix.NatomE(self.E.ind_O); end
if isempty(self.E.ind_N), mix.w = 0; else, mix.w = mix.NatomE(self.E.ind_N); end

mix.N      = sum(SpeciesMatrix(:, 1)); %[mol] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MARCOS: rho and mass fractions computations
% Ni = SpeciesMatrix(:,1); % number of moles of species i in the mixture
% 
% NM = sum(Ni);   % overall number of moles in the mixture
% Xi = Ni/NM;     % mass fraction of species i
% 
% MM = 0;         % overall mass of mixture
% WM = 0;         % averge moleular mass of mixture
% 
% Yi_times_WM = zeros(size(Xi));
% 
% fnm = fieldnames(DB);
% for i = 1:numel(fnm)
%     Wi = DB.(fnm{i}).mm/1000;
%     MM = MM + Wi*Ni(i);
%     WM = WM + Wi*Xi(i);
%     Yi_times_WM(i) = Wi*Xi(i);
% end
% str.MM = MM;
% Yi = Yi_times_WM/WM;
% Rg = R0/WM;
% str.W2 = WM*1e3;
% str.rho2 = p*101325/(Rg*T);
% str.Xi = Xi';
% str.Yi = Yi';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% str.hf     = sum(SpeciesMatrix(:,2));
% str.DhT    = sum(SpeciesMatrix(:,3));
% str.h_marcos      = str.hf + str.DhT;
% str.hf     = sum(SpeciesMatrix(:,2)./Ni.*Yi,'OmitNan');
% str.DhT     = sum(SpeciesMatrix(:,3)./Ni.*Yi,'OmitNan');
% str.h      = str.hf + str.DhT;
% str.ef     = sum(SpeciesMatrix(:,4)./Ni.*Yi,'OmitNan');
% str.DeT    = sum(SpeciesMatrix(:,5)./Ni.*Yi,'OmitNan');
% str.e      = str.ef + str.DeT;
% str.cP     = sum(SpeciesMatrix(:,6)./Ni.*Yi,'OmitNan');    
% str.cV     = sum(SpeciesMatrix(:,7)./Ni.*Yi,'OmitNan');    
% str.pv     = sum(SpeciesMatrix(:,9));
% str.p      = p;
% str.v      = str.pv/str.p;
% swtCond    = SpeciesMatrix(:,10);
% str.mi     = sum(SpeciesMatrix(:,11));
% str.rho    = str.mi/str.v;
% 
% Yi         = SpeciesMatrix(:,11)./str.mi;
% str.W      = 1/sum(Yi./SpeciesMatrix(:,12),'OmitNan');
% 
% Ngas       = sum(SpeciesMatrix(:,1).*(1-swtCond));
% 
% Ni         = SpeciesMatrix(:,1);
% Xi         = Ni/Ngas;
% ii         = find(Xi>0);
% 
% str.S0     = sum(SpeciesMatrix(:,8)./Ni.*Yi,'OmitNan');
% DSi        = -R0*Ni(ii).*log(Xi(ii)*p).*(1-swtCond(ii)); % only nonzero for noncondensed species
% str.DS     = sum(DSi);
% str.S      = str.S0+str.DS;


%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
mix.hf     = sum(SpeciesMatrix(:,2));  % [kJ]
mix.DhT    = sum(SpeciesMatrix(:,3));  % [kJ]
mix.h      = mix.hf + mix.DhT;         % [kJ]
mix.ef     = sum(SpeciesMatrix(:,4));  % [kJ]
mix.DeT    = sum(SpeciesMatrix(:,5));  % [kJ]
mix.e      = mix.ef + mix.DeT;         % [kJ]
mix.cP     = sum(SpeciesMatrix(:,6));  % [J/K]     
mix.cV     = sum(SpeciesMatrix(:,7));  % [J/K]   
mix.pv     = sum(SpeciesMatrix(:,9));
mix.p      = p;
mix.v      = mix.pv/mix.p;
mix.swtCond = SpeciesMatrix(:,10);
mix.mi     = sum(SpeciesMatrix(:,11))*1e-3; % [kg]
mix.rho    = mix.mi/mix.v*1e3;              % [kg/m3]

mix.Yi     = SpeciesMatrix(:,11)./mix.mi*1e-3;   % [-]
mix.W      = 1/sum(mix.Yi./SpeciesMatrix(:,12),'OmitNan');

% Ngas       = sum(SpeciesMatrix(:,1));
% Ngas       = sum(SpeciesMatrix(:,1).*(1-str.swtCond));
% htest = sum(str.Yi .* str.hf) + sum(str.Yi .* str.Dht)
Ni         = SpeciesMatrix(:,1); % [mol]
% str.Xi     = Ni/Ngas;               % [-]
mix.Xi     = Ni/mix.N;               % [-]
ii         = mix.Xi>0;     
mix.S0     = sum(SpeciesMatrix(:,8)); % [kJ/K]
DSi        = Ni(ii).*log(mix.Xi(ii)*p).*(1-mix.swtCond(ii)); % only nonzero for noncondensed species
mix.DS     = -R0*sum(DSi) * 1e-3; % [kJ/K]
mix.S      = (mix.S0 + mix.DS); % [kJ/K]
mix.T = T; % [K]
mix.g = mix.h - mix.T * mix.S; % [kJ]
mix.e = mix.h - sum(Ni(ii).*(1-mix.swtCond(ii)))*R0*T*1e-3; % THIS WORKS!!!
% str.e = str.h - sum(Ni(ii))*R0*T*1e-3; % THIS WORKS!!!
mix.gamma = mix.cP/mix.cV;
mix.sound = sqrt(mix.gamma*p*1e5/mix.rho);
% Correction of: cP, cV, gamma and speed of sound as consequence of the
% chemical reaction
if isfield(self, 'dNi_T')
    ind = SpeciesMatrix(:, 1) > 0;
    H0_j = (SpeciesMatrix(ind, 2) + SpeciesMatrix(ind, 3)) * 1e3;
    mix.cP_r = sum(H0_j/T .* self.dNi_T(ind));
    mix.cP = mix.cP + mix.cP_r;
    mix.dVdT_p =  1 + self.dN_T;
    mix.dVdp_T = -1 + self.dN_P;
    mix.cV = mix.cP + (mix.pv/T * mix.dVdT_p^2) / mix.dVdp_T *1e2;
    mix.gamma = mix.cP/mix.cV;
    mix.gamma_s = - mix.gamma / mix.dVdp_T;
    mix.sound = sqrt(mix.gamma_s*p*1e5/mix.rho);
end


