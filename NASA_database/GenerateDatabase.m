function strThProp = GenerateDatabase(strMaster)
if ~exist('strThProp', 'var')
    if exist('strThProp.mat', 'file')
        fprintf('NASA short database loaded from main path ... ')
        load('strThProp.mat', 'strThProp');
    else
        strThProp = get_Database(strMaster); % struct with tabulated data of selected species.
    end
    fprintf('OK!\n');
else
    fprintf('NASA short database already loaded\n')
end
end
function strThProp = get_Database(strMaster)

%%% gri30-x.cti except 'CH2(s)' + others
SpeciesList = {'H2', 'H', 'O', 'O2', 'O3', 'OH', 'H2O', 'HO2', 'H2O2', 'C',...
    'CH', 'CH2', 'CH3', 'CH4', 'CO', 'CO2', 'HCO',...
    'CH2OH', 'CH3O', 'CH3OH', 'C2H', 'C2H4',...
    'C2H5', 'C2H6', 'HCCO', 'N', 'NH', 'NH2', 'NH3', 'N2H4',...
    'NO', 'NO2', 'NO3', 'N2O', 'N2O3', 'N2O4', 'HNO', 'CN', 'HCN',...
    'NCO', 'N2', 'N3', 'Ar', 'C3H8','C2','C2H2_acetylene','C6H6',...
    'C8H18_isooctane','C2H5OH','He','Cbgrb'};

% SpeciesList = {'O2','N2','O','O3','N','NO','NO2','NO3','N2O','N2O3','N2O4','N3', ...
%     'eminus', 'Arplus', 'Nminus', 'Nplus', 'NOplus', 'NO2minus', 'NO3minus', 'N2plus', 'N2minus', 'N2Oplus', ...
%      'Oplus', 'Ominus', 'O2plus', 'O2minus'};
 

% SpeciesList = {'O2','N2','O','N'};

% SpeciesList = {'N2', 'N', 'N3', 'eminus', 'N2plus', 'N2minus', 'Nplus', 'Nplus'};

% SpeciesList = {'N2', 'N', 'N3'};

% SpeciesList = {'O3', 'O2', 'O', 'Oplus', 'Ominus', 'O2plus', 'O2minus', 'eminus'};

% SpeciesList = {'H2', 'H', 'Hplus', 'Hminus', 'H2plus', 'H2minus', 'eminus'};

% SpeciesList = {'CO2', 'CO', 'H2O', 'H2', 'O2', 'N2', 'He', 'Ar',...
%                 'OH','H','O','HO2','NO','HCO','CH4','CH3',...
%                 'NO2','NH3','NH2','N','HCN','CN','N2O','C2','CH',...
%                 'H2bLb','O2bLb','RP_1'};

%  NASA *: CH4 + 2O2 + 7.52N2
% SpeciesList = {'C','CN','CO2','H','O2','CH','C3','C5','NH','O','CO',...
%     'C2','C4','H2','N','NO','N2','OH','CH4','H2O','He','Ar','Cbgrb'};


% NASA ALL: CH4 + 2O2 + 7.52N2
% EXCEPTIONS: 'THDCPD_exo', 'C6H14bLb_n_hexa','H2Obcrb','bHCOOHb2' ,'bCH3COOHb2'
% SpeciesList = {'C','CH3','CH4','CN','CO2','C2H','CH2CO_ketene',...
%     'C2H3_vinyl','C2H4','CH3COOH','C2H6','CH3OCH3','CNC','C2O',...
%     'C3H3_2_propynl','C3H6O_propanal',...
%     'C3H8','CNCOCN','C4H2_butadiyne','C4H6_1butyne','C4H8_1_butene',...
%     'C4H8_isobutene','C4H9_n_butyl','C4H9_t_butyl','C4N2',...
%     'C5H11_pentyl','C5H12_i_pentane','C6H5_phenyl','C6H5OH_phenol',...
%     'C7H7_benzyl','C7H14_1_heptene','C7H16_2_methylh',...
%     'C8H16_1_octene','C8H18_isooctane','C10H21_n_decyl','H','HCCN','HNCO',...
%     'HNO3','HCHO_formaldehy','H2O2','NCO','NH3','NO2','N2O','NH2NO2','N2O4',...
%     'N3H','O2','Cbgrb','C2H5OHbLb','C6H6bLb','H2ObLb','CH','CH2OH','CH3OH',...
%     'CNN','COOH','C2H2_acetylene','ObCHb2O','CH3CN','C2H4O_ethylen_o',...
%     'OHCH2COOH','CH3N2CH3','CH3O2CH3','OCCN','C3','C3H4_allene','C3H5_allyl',...
%     'C3H6O_propylox','C3H7_n_propyl','C3H8O_1propanol','C3O2',...
%     'C4H6_2butyne','C4H8_cis2_buten',...
%     'C4H9_i_butyl','C4H10_n_butane','C5','C5H10_1_pentene','C5H11_t_pentyl',...
%     'CH3CbCH3b2CH3','C6H5O_phenoxy','C6H13_n_hexyl','C7H8',...
%     'C7H15_n_heptyl','C8H8_styrene','C8H17_n_octyl','C9H19_n_nonyl','C12H9_o_bipheny',...
%     'HCN','HCCO','HNO','HO2','HCOOH','NH','NH2OH','NO3','NCN','N2H4',...
%     'N2O5','O','O3','N2H4bLb','CH2',...
%     'CH3O','CH3OOH','CO','C2','C2H2_vinylidene','HObCOb2OH','CH3CO_acetyl',...
%     'CH3CHO_ethanal','C2H5','C2H5OH','CCN','C2N2','C3H3_1_propynl','C3H4_propyne',...
%     'C3H6_propylene','C3H6O_acetone','C3H7_i_propyl','C3H8O_2propanol',...
%     'C4','C4H6_butadiene','C4H8_tr2_butene',...
%     'C4H9_s_butyl','C4H10_isobutane',...
%     'C5H12_n_pentane','C6H2','C6H6','C6H12_1_hexene','C6H14_n_hexane',...
%     'C7H8O_cresol_mx','C7H16_n_heptane','C8H10_ethylbenz','C8H18_n_octane',...
%     'C10H8_naphthale','C12H10_biphenyl','HCO','HNC','HNO2','H2','H2O',...
%     'N','NH2','NO','N2','N2H2','N2O3','N3','OH','CH3OHbLb','C6H5NH2bLb',...
%     'He','Ar','Cbgrb'};

% SHOCK NASA O2+N2 + OTHERS
% SpeciesList = {'O2','N2','O','O3','N','NO','NO2','NO3','N2O','N2O3','N2O4','N3',...
%     'C','C2','CO','CO2','CN','Ar','CH4','H2O','H2','H','He','OH','Cbgrb'};

% HYDROGEN
% SpeciesList = {'H','HNO','HNO3','H2O','NH','NH2OH','NO3','N2H2','N2O3','N3','OH','HNO2',...
%                 'H2','N','NH3','NO2','N2O','N2H4','N2O5','O','O3','O2','N2','HO2','NH2','H2O2',...
%                 'N3H','NH2NO2'};

% SpeciesList = {'RP_1'};

% SpeciesList = {'O2','N2','O','O3','N','NO','NO2','NO3','N2O','N2O3','N2O4','N3'};

% SpeciesList = {'O2','N2','O','O3','N','NO','NO2','NO3','N2O','N2O3','N2O4','N3',...
%     'C','C2','CO','CO2','CN','Ar','CH4','H2O','H2','H','He','OH','Cbgrb','F','F2'};

fprintf('Generating short NASA database ... ')
for i = 1:length(SpeciesList)
    
    Species = FullName2name(SpeciesList{i});
    
    if isfield(strMaster,Species)
        
        % disp(SpeciesList{i})
        
        ctTInt = strMaster.(Species).ctTInt;
        tRange = strMaster.(Species).tRange;
        swtCondensed = sign(strMaster.(Species).swtCondensed);
        
        if ctTInt > 0
            
            [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = SpeciesThermProp(strMaster,SpeciesList{i},298.15,'molar',0);
            
            strThProp.(Species).name = Species;
            strThProp.(Species).FullName = SpeciesList{i};
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
                [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = SpeciesThermProp(strMaster,SpeciesList{i},T,'molar',0);
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
            strThProp.(Species).cPcurve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).cp, 'pchip', 'linear');
            strThProp.(Species).cVcurve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).cv, 'pchip', 'linear');
            strThProp.(Species).DhTcurve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).DhT, 'pchip', 'linear');
            strThProp.(Species).DeTcurve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).DeT, 'pchip', 'linear');
            strThProp.(Species).h0curve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).h0, 'pchip', 'linear');
            strThProp.(Species).s0curve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).s0, 'pchip', 'linear');
            strThProp.(Species).g0curve = griddedInterpolant(strThProp.(Species).T,strThProp.(Species).g0, 'pchip', 'linear');
            
            % DATA COEFFICIENTS NASA 9 POLYNOMIAL (ONLY GASES)
            strThProp.(Species).ctTInt = strMaster.(Species).ctTInt;
            strThProp.(Species).tRange = strMaster.(Species).tRange;
            strThProp.(Species).tExponents = strMaster.(Species).tExponents;
            strThProp.(Species).ctTInt = strMaster.(Species).ctTInt;
            strThProp.(Species).a = strMaster.(Species).a;
            strThProp.(Species).b  = strMaster.(Species).b;
        else
            
            Tref = tRange(1);
            
            [txFormula, mm, Cp0, Cv0, Hf0, H0, Ef0, E0, S0, DfG0] = SpeciesThermProp(strMaster,SpeciesList{i},Tref,'molar',0);
            
            strThProp.(Species).name = Species;
            strThProp.(Species).FullName = SpeciesList{i};
            strThProp.(Species).txFormula = txFormula;
            strThProp.(Species).mm  = mm;
            strThProp.(Species).hf  = Hf0;
            strThProp.(Species).ef  = Hf0;
            strThProp.(Species).swtCondensed = swtCondensed;
            strThProp.(Species).T   = Tref;
            strThProp.(Species).DhT = 0;
            strThProp.(Species).DeT = 0;
            strThProp.(Species).h0  = H0;
            strThProp.(Species).s0  = S0;
            strThProp.(Species).cp  = Cp0;
            strThProp.(Species).cv  = Cv0;
            strThProp.(Species).g0  = DfG0;
            
            strThProp.(Species).ctTInt = 0;
        end
    else
        fprintf(['\n- Species ''', SpeciesList{i}, ''' does not exist as a field in strMaster structure ... '])
    end
    
end
end