function mix2 = equilibrateT(obj, mix1, mix2, T, varargin)
    % Obtain equilibrium properties and composition for the given
    % temperature [K] and pressure [bar] / specific volume [m3/kg]
    %
    % Args:
    %     obj (EquilibriumSolver): EquilibriumSolver object
    %     mix1 (Mixture): Properties of the initial mixture
    %     mix2 (Mixture): Properties of the final mixture
    %     T (float): Temperature [K]
    %
    % Optional Args:
    %     * molesGuess (float): Mixture composition [mol] of a previous computation
    %
    % Returns:
    %     mix2 (Mixture): Properties of the final mixture
    %
    % Examples:
    %     mix2 = equilibrateT(EquilibriumSolver(), mix1, mix2, 3000)
    %     mix2 = equilibrateT(EquilibriumSolver(), mix1, mix2, 3000, molesGuess)
    
    % Check if calculations are for a calorically perfect gas
    if obj.caloricGasModel.isPerfect()
        mix2 = equilibrateTPerfect(mix1, mix2, T);
        return
    end

    % Check if calculations are for a thermally perfect gas
    if obj.caloricGasModel.isThermallyPerfect()
        mix2 = equilibrateTFrozen(mix2, T);
        return
    end

    % Definitions
    molesGuess = [];
    N_mix0 = moles(mix1); % Get moles of inert species
    system = mix2.chemicalSystem;
    productSpeciesSet = mix2.productSpeciesSet;
    
    % Unpack addtitional inputs
    if nargin > 4
        molesGuess = varargin{1};
        
        if system.FLAG_COMPLETE
            molesGuess = [];
        elseif ~isempty(molesGuess)
            molesGuess = molesGuess(productSpeciesSet.indexGlobal);
        end

    end

    % Check flag
    if ~obj.FLAG_FAST, molesGuess = []; end

    % Compute number of moles
    [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = selectEquilibrium(obj, system, productSpeciesSet, T, mix1, mix2, molesGuess);
    
    % Compute mixture state at chemical equilibrium
    setMixture(mix2);

    % NESTED FUNCTIONS
    function mix = setMixture(mix)
        % Compute properties of final mixture

        % Scatter only the final active equilibrium species; candidate product slots may contain trial values
        indexActive = productSpeciesSet.indexGlobal(indexProducts);
        N = productVectorToSystemVector(system.numSpecies, productSpeciesSet, N, indexProducts);
        dNi_T = productVectorToSystemVector(system.numSpecies, productSpeciesSet, dNi_T, indexProducts);
        dNi_p = productVectorToSystemVector(system.numSpecies, productSpeciesSet, dNi_p, indexProducts);
        h0 = productVectorToSystemVector(system.numSpecies, productSpeciesSet, h0, indexProducts);
        N(system.indexFrozen) = N_mix0(system.indexFrozen);

        % Assign values
        mix.T = T;
        mix.dNi_T = dNi_T; mix.dN_T = dN_T;
        mix.dNi_p = dNi_p; mix.dN_p = dN_p;
        mix.FLAG_REACTION = true;
        mix.errorMoles = STOP;
        mix.errorMolesIons = STOP_ions;

        % Compute properties of final mixture
        setMolesFast(mix, N, h0, indexActive);
    end
    
end

% SUB-PASS FUNCTIONS
function [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = selectEquilibrium(obj, system, productSpeciesSet, T, mix1, mix2, molesGuess)
    % Select equilibrium: TP: Gibbs; TV: Helmholtz
    if strfind(obj.problemType, 'P') == 2
        [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = equilibriumGibbs(obj, system, productSpeciesSet, mix2.p, T, mix1, molesGuess);
        return
    end
    
    [N, dNi_T, dN_T, dNi_p, dN_p, indexProducts, STOP, STOP_ions, h0] = equilibriumHelmholtz(obj, system, productSpeciesSet, mix2.v, T, mix1, molesGuess);
end

function vector = productVectorToSystemVector(numSpecies, productSpeciesSet, vectorProduct, indexProduct)
    % Map product species vector into ChemicalSystem species order
    vector = zeros(numSpecies, 1);
    vector(productSpeciesSet.indexGlobal(indexProduct)) = vectorProduct(indexProduct);
end

function mix2 = equilibrateTPerfect(mix1, mix2, T)
    % Obtain equilibrium properties and composition for the given
    % temperature [K] and pressure [bar] assuming a calorically perfect gas
    %
    % Args:
    %     mix1 (Mixture): Properties of the initial mixture
    %     mix2 (Mixture): Properties of the final mixture
    %     T (float): Temperature [K]
    %
    % Returns:
    %     mix2 (Mixture): Properties of the final mixture assuming a calorically perfect gas
    %
    % Example:
    %     mix2 = equilibrateTPerfect(mix1, mix2, 3000)
    
    % Import packages
    import combustiontoolbox.common.Units

    % Definitions
    Tref = 298.15; % [K] To be fixed: each species should have its own reference temperature

    % Recompute properties of mix2
    setTemperature(mix2, T);

    % Change properties that remains thermochemically frozen
    mix2.cp = mix1.cp;
    mix2.cv = mix1.cv;
    mix2.gamma = mix1.gamma;
    mix2.gamma_s = mix1.gamma_s;
    mix2.sound = sqrt(mix2.gamma * mix2.p * Units.bar2Pa / mix2.rho);

    % Compute enthalpy [J]
    mix2.hf = mix1.hf;
    mix2.DhT = mix1.cp * (T - Tref);
    mix2.h = mix2.hf + mix2.DhT;

    % Compute internal energy [J]
    mix2.ef = mix1.ef;
    mix2.DeT = mix1.cv * (T - Tref);
    mix2.e = mix2.ef + mix2.DeT;

    % Compute entropy [J/K]
    mix2.s0 = mix1.s0 + mix1.cp * log(T / mix1.T);
    mix2.s = mix2.s0 + mix2.Ds;
    
    % Compute Gibbs free energy
    mix2.g = mix2.h - T * mix2.s;

    if isempty(mix2.u)
        return
    end
    
    mix2.uShock = mix2.u;
    mix2.mach = mix2.u / mix2.sound;
end

function mix = equilibrateTFrozen(mix, T)
    % Obtain equilibrium properties and composition for the given
    % temperature [K] and pressure [bar] assuming a thermally perfect gas
    %
    % Args:
    %     mix (Mixture): Properties the mixture
    %     T (float): Temperature [K]
    %
    % Returns:
    %     mix (Mixture): Properties of the final mixture assuming a thermally perfect gas
    %
    % Example:
    %     mix = equilibrateTFrozen(mix, 3000)

    % Definitions
    mix.FLAG_REACTION = false;
    mix.T = T;
    mix.quantity = mix.Xi * mix.N;
    mix.listSpecies = mix.chemicalSystem.listSpecies;
    
    % Update indexSpecies
    updateIndexSpecies(mix);

    % Update thermodynamic properties assuming a thermally perfect gas
    updateThermodynamics(mix); 
end