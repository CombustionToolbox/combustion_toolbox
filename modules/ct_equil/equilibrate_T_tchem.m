function mix2 = equilibrate_T_tchem(self, mix1, pP, TP)
    % Obtain equilibrium properties and composition for the given
    % temperature [K] and pressure [bar] assuming a calorically frozen gas
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     mix1 (struct): Properties of the initial mixture
    %     pP (float): Pressure [bar]
    %     TP (float): Temperature [K]
    %
    % Returns:
    %     mix2 (struct): Properties of the final mixture assuming a calorically frozen gas
    %
    % Example:
    %     mix2 = equilibrate_T_tchem(self, self.PS.strR{1}, 1.01325, 3000)

    % Initialize mix2 as mix1 (most properties are the same)
    mix2 = mix1;
    % Set temperature [K] of the mixture 
    mix2.T = TP;
    % Set thermodynamic derivative to their frozen values
    mix2.dVdT_p = 1;
    mix2.dVdp_T = -1;
    % Set pressure and volume of the mixture
    if strfind(self.PD.ProblemType, 'P') == 2
        % Set pressure [bar] of the mixture
        mix2.p = pP;
        % Set volume [m3] of the mixture
        mix2.v = self.PD.EOS.volume(self, TP, convert_bar_to_Pa(pP), self.S.LS, mix2.Xi) * molesGas(mix2);
    else
        % Set volume [m3] of the mixture
        mix2.v = mix1.v; % Included for clarification
        % Set pressure [bar] of the mixture
        mix2.p = convert_Pa_to_bar(self.PD.EOS.pressure(self, molesGas(mix2), TP, mix2.v));
    end
    % Set density [kg/m3] of the mixture
    mix2.rho = mix2.mi / mix2.v;
    % Set sound velocity [m/s]
    mix2.sound = sqrt(mix2.gamma * convert_bar_to_Pa(mix2.p) / mix2.rho);
end