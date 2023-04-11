function HTR_shock_incident(T1, p1, M1)
    % Get input conditions for the Hypersonic Task-based Research (HTR) solver [1]
    % (see the GitHub repository [2]) using the Combustion Toolbox [3].
    %
    % References:
    %     [1] Di Renzo, M., Fu, L., and Urzay, J., HTR solver: An open-source
    %         exascale-oriented task-based multi-GPU high-order code for
    %         hypersonics aerothermodynamics, Comput. Phys. Commun, Vol. 255,
    %         2020, p. 107262
    %
    %     [2] https://github.com/stanfordhpccenter/HTR-solver
    %
    %     [3] Cuadra, A., Huete, C., Vera, M., (2022). Combustion Toolbox:
    %         A MATLAB-GUI based open-source tool for solving gaseous
    %         combustion problems. Zenodo. DOI: 10.5281/zenodo.5554911.
    %
    % @author: Alberto Cuadra Lara
    %          PhD Candidate - Group Fluid Mechanics
    %          Universidad Carlos III de Madrid
    %                  
    % Last update Apr 05 2023

    % Inputs
    Oxidizer = {'N2', 'O2'};
    moles = [79, 21] / 21;
    LS =  {'N2', 'O2'};
    
    % Initialize CT
    self = App();

    % Compute frozen sound velocity [m/s]
    a1 = compute_sound(T1, p1, Oxidizer, moles, 'self', self);

    % Compute pre-shock velocity [m/s]
    u1 = a1 * M1;

    % Tuning parameters
    tolN = 1e-14;

    % Combustion Toolbox
    results_CT = run_CT('ProblemType', 'SHOCK_I',...
                        'Species', LS,...
                        'S_Oxidizer', Oxidizer,...
                        'N_Oxidizer', moles,...
                        'TR', T1,...
                        'PR', p1,...
                        'u1', u1,...
                        'tolN', tolN,...
                        'DB', self.DB,...
                        'DB_master', self.DB_master);

    % Combustion Toolbox - thermochemically frozen (calorically perfect gas)
    results_CT_tchem = run_CT('ProblemType', 'SHOCK_I',...
                               'Species', Oxidizer,...
                               'S_Oxidizer', Oxidizer,...
                               'N_Oxidizer', moles,...
                               'TR', T1,...
                               'PR', p1,...
                               'u1', u1,...
                               'tolN', tolN,...
                               'FLAG_TCHEM_FROZEN', true,...
                               'DB', self.DB,...
                               'DB_master', self.DB_master);
    
    % Get and print results
    get_results(results_CT, u1, false);
    get_results(results_CT_tchem, u1, true);
end

% SUB-PASS FUNCTIONS
function get_results(self, u1, FLAG_FROZEN)
    % Get results from self

    % Get list and number of species
    LS = self.S.LS;
    NS = self.S.NS;

    % Get mixtures
    mix1 = self.PS.strR;
    mix2 = self.PS.strP;

    % Define gas constant
    Rgas = self.C.R0 / (MolecularWeight(mix1{1}) * 1e-3);

    % Get results
    rho1 = cell2vector(mix1, 'rho');
    rho2 = cell2vector(mix2, 'rho');
    T1 = cell2vector(mix1, 'T');
    T2 = cell2vector(mix2, 'T');
    p1 = cell2vector(mix1, 'p');
    p2 = cell2vector(mix2, 'p');
    sound1 = cell2vector(mix1, 'sound');
    sound2 = cell2vector(mix2, 'sound');
    v_shock = cell2vector(mix2, 'v_shock');
    Xi = cell2vector(mix2, 'Xi');
    
    % Set reference values
    p_ref = convert_bar_to_Pa(p1); % [Pa]
    T_ref = T1; % [K]
    rho_ref = p_ref / (Rgas * T_ref);
    u_ref = sqrt(p_ref / rho_ref);

    % Compute Mach numbers
    M1 = u1 ./ sound1;
    M2 = v_shock ./ sound2;

    % Compute ratios
    R = rho2 / rho1;
    T = T2 / T1;
    P = p2 / p1;
    U = M1 * sound1 / u_ref;
    
    % Display results
    fprintf('\n********************************************\n');
    if FLAG_FROZEN
        fprintf('CASE: THERMOCHEMICALLY FROZEN\n');
    else
        fprintf('CASE: FROZEN\n');
    end
    
    if ~FLAG_FROZEN
        fprintf('Case for HTR solver:\n\n');
        fprintf('Pre-shock Mach       |  %5.4f\n', M1);
        fprintf('Post-shock vel [m/s] |  %5.4f\n', v_shock);
        fprintf('Pre-shock sound      |  %5.4f\n', sound1);
        fprintf('Post-shock sound     |  %5.4f\n\n', sound2);
    end
    
    fprintf('Initialization values for HTR solver:\n\n');
    fprintf('Pressure ratio       |  %5.4f\n', P);
    fprintf('Temperature ratio    |  %5.4f\n', T);
    fprintf('Downstream Mach      |  %5.4f\n', M2);
    fprintf('Velocity ratio       |  %5.4f\n', U);
    fprintf('Molar fractions      |\n');
    fprintf('   %16s  |  %1.4e\n', LS{1}, Xi(1));
    fprintf('   %16s  |  %1.4e\n', LS{2}, Xi(2));
    if NS > 2
        fprintf('   %16s  |  %1.4e\n', LS{3}, Xi(3));
        fprintf('   %16s  |  %1.4e\n', LS{4}, Xi(4));
        fprintf('   %16s  |  %1.4e\n', LS{5}, Xi(5));
        fprintf('\nXi : [%1.4e, %1.4e, %1.4e, %1.4e, %1.4e]\n', Xi);
    end
end