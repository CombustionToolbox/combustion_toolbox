function str =  ComputeProperties(A0,SpeciesMatrix,p,T,ind_C,ind_H)
R0 = 8.3144598; % [J/(K mol)]. Universal gas constant
str.NatomE = sum(SpeciesMatrix(:,1).*A0);
str.x     = str.NatomE(ind_C); % NatomE(find(strcmp(Elements,'C')));
str.y     = str.NatomE(ind_H); % NatomE(find(strcmp(Elements,'H')));

str.N      = sum(SpeciesMatrix(:,1)); %[mol]
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
% fnm = fieldnames(strThProp);
% for i = 1:numel(fnm)
%     Wi = strThProp.(fnm{i}).mm/1000;
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
str.hf     = sum(SpeciesMatrix(:,2));  % [kJ]
str.DhT    = sum(SpeciesMatrix(:,3));  % [kJ]
str.h      = str.hf + str.DhT;         % [kJ]
str.ef     = sum(SpeciesMatrix(:,4));  % [kJ]
str.DeT    = sum(SpeciesMatrix(:,5));  % [kJ]
str.e      = str.ef + str.DeT;         % [kJ]
str.cP     = sum(SpeciesMatrix(:,6));  % [J/K]    
str.cV     = sum(SpeciesMatrix(:,7));  % [J/K]   
str.pv     = sum(SpeciesMatrix(:,9));
str.p      = p;
str.v      = str.pv/str.p;
str.swtCond = SpeciesMatrix(:,10);
str.mi     = sum(SpeciesMatrix(:,11))*1e-3; % [kg]
str.rho    = str.mi/str.v*1e3;                 % [kg/m^3]

str.Yi     = SpeciesMatrix(:,11)./str.mi*1e-3;   % [-]
str.W      = 1/sum(str.Yi./SpeciesMatrix(:,12),'OmitNan');

% Ngas       = sum(SpeciesMatrix(:,1));
% Ngas       = sum(SpeciesMatrix(:,1).*(1-str.swtCond));

Ni         = SpeciesMatrix(:,1); % [mol]
% str.Xi     = Ni/Ngas;               % [-]
str.Xi     = Ni/str.N;               % [-]
ii         = find(str.Xi>0);         

str.S0     = sum(SpeciesMatrix(:,8)); % [kJ/K]
% DSi        = -R0*Ni(ii).*log(str.Xi(ii)*p).*(1-swtCond(ii)); % only nonzero for noncondensed species
DSi        = Ni(ii).*log(str.Xi(ii)*p).*(1-str.swtCond(ii)); % only nonzero for noncondensed species
str.DS     = -R0*sum(DSi); % [J/K]
str.S      = (str.S0+str.DS*1e-3)/str.mi; % [kJ/(kg-K)]
str.T = T; % [K]
% str.e = str.h - sum(Ni(ii).*(1-str.swtCond(ii)))*R0*T*1e-3; % THIS WORKS!!!
str.e = str.h - sum(Ni(ii))*R0*T*1e-3; % THIS WORKS!!!
str.gamma = str.cP/str.cV;
str.sound = sqrt(str.gamma*p*1e5/str.rho);