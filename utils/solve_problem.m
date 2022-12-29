function self = solve_problem(self, ProblemType)
    % Solve the given ProblemType with the conditions and mixture specified in self
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %     ProblemType (str): Tag of the problem to solve
    %
    % Returns:
    %     self (struct): Data of the mixtures (initial and final), conditions, databases

    try
        % Set Problem Type
        self.PD.ProblemType = ProblemType;
        % Check inputs and set length of the loop
        self = check_inputs(self);
        % Get Flags and length of the loop
        self = get_FLAG_N(self);
        % Loop
        for i = self.C.l_phi:-1:1
            % Set problem conditions by case
            self = set_problem_conditions(self, i);
            % Compute the properties of the initial mixture
            self = define_FOI(self, i);
            % Check if list of products corresponds to "complete reaction"
            self = check_complete_reaction(self, i);
            % Solve selected problem
            self = select_problem(self, i);
            % Display results in the command window
            results(self, i);
        end

        % Print execution time
        self.Misc.timer_loop = toc(self.Misc.timer_0);
        fprintf('Elapsed time is %.6g seconds\n', self.Misc.timer_loop);
        % Retrieve range values (problem set only)
        if isfield(self.PD, 'range_name')
            self.PD.(self.PD.range_name).value = self.PD.range;
        end

    catch ME
        print_error(ME);
    end

end

% SUB-PASS FUNCTIONS
function self = set_problem_conditions(self, i)
    % Set problem conditions per case
    if isfield(self.PD, 'range_name')

        if ~strcmpi(self.PD.range_name, 'phi')
            self.PD.(self.PD.range_name).value = self.PD.range(i);
        end

    end

end

function self = select_problem(self, i)
    % Solve selected problem
    switch upper(self.PD.ProblemType)
        case {'SHOCK_I', 'SHOCK_R'}

            try
                u1 = self.PD.u1.value(i);
            catch
                u1 = self.PD.u1.value;
            end

            if strcmp(self.PD.ProblemType, 'SHOCK_I')

                if i == self.C.l_phi
                    [self.PS.strR{i}, self.PS.strP{i}] = shock_incident(self, self.PS.strR{i}, u1);
                else
                    [self.PS.strR{i}, self.PS.strP{i}] = shock_incident(self, self.PS.strR{i}, u1, self.PS.strP{i + 1});
                end

            else

                if i == self.C.l_phi
                    % Calculate post-shock state (2)
                    [self.PS.strR{i}, self.PS.str2{i}] = shock_incident(self, self.PS.strR{i}, u1);
                    % Calculate post-shock state (5)
                    [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, u1, self.PS.str2{i});
                else
                    % Calculate post-shock state (2)
                    [self.PS.strR{i}, self.PS.str2{i}] = shock_incident(self, self.PS.strR{i}, u1, self.PS.str2{i + 1});
                    % Calculate post-shock state (5)
                    [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, u1, self.PS.str2{i}, self.PS.strP{i + 1});
                end

            end

        case {'SHOCK_OBLIQUE'}

            try
                u1 = self.PD.u1.value(i);
            catch
                u1 = self.PD.u1.value;
            end

            if ~isempty(self.PD.theta.value)

                try
                    theta = self.PD.theta.value(i);
                catch
                    theta = self.PD.theta.value;
                end

                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_oblique_theta(self, self.PS.strR{i}, u1, theta);
            else

                try
                    beta = self.PD.beta.value(i);
                catch
                    beta = self.PD.beta.value;
                end

                [self.PS.strR{i}, self.PS.strP{i}] = shock_oblique_beta(self, self.PS.strR{i}, u1, beta);
            end

        case {'SHOCK_OBLIQUE_R'}

            try
                u1 = self.PD.u1.value(i);
            catch
                u1 = self.PD.u1.value;
            end

            if ~isempty(self.PD.theta.value)

                try
                    theta = self.PD.theta.value(i);
                catch
                    theta = self.PD.theta.value;
                end

                self.PD.ProblemType = 'SHOCK_OBLIQUE';
                [self.PS.strR{i}, self.PS.str2{i}, ~] = shock_oblique_theta(self, self.PS.strR{i}, u1, theta);
                self.PD.ProblemType = 'SHOCK_OBLIQUE_R';
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.str3{i}, self.PS.str3_2{i}] = shock_oblique_reflected_theta(self, self.PS.strR{i}, self.PS.str2{i}.u, self.PS.str2{i}.theta, self.PS.str2{i});
            else

                try
                    beta = self.PD.beta.value(i);
                catch
                    beta = self.PD.beta.value;
                end

                self.PD.ProblemType = 'SHOCK_OBLIQUE';
                [self.PS.strR{i}, self.PS.str2{i}] = shock_oblique_beta(self, self.PS.strR{i}, u1, beta);
                self.PD.ProblemType = 'SHOCK_OBLIQUE_R';
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.str3{i}, self.PS.str3_2{i}] = shock_oblique_reflected_theta(self, self.PS.strR{i}, self.PS.str2{i}.u, self.PS.str2{i}.theta, self.PS.str2{i});
            end

        case {'SHOCK_POLAR'}

            try
                u1 = self.PD.u1.value(i);
            catch
                u1 = self.PD.u1.value;
            end

            [self.PS.strR{i}, self.PS.strP{i}] = shock_polar(self, self.PS.strR{i}, u1);
        case {'SHOCK_POLAR_R'}

            try
                u1 = self.PD.u1.value(i);
            catch
                u1 = self.PD.u1.value;
            end

            if ~isempty(self.PD.theta.value)

                try
                    theta = self.PD.theta.value(i);
                catch
                    theta = self.PD.theta.value;
                end

                [self.PS.strR{i}, self.PS.str2{i}] = shock_polar(self, self.PS.strR{i}, u1);
                [~, self.PS.str2_1{i}] = shock_oblique_theta(self, self.PS.strR{i}, u1, theta);
                [~, self.PS.strP{i}] = shock_polar(self, self.PS.str2_1{i}, self.PS.str2_1{i}.u);
                [~, self.PS.str3_1{i}, self.PS.str3_2{i}] = shock_oblique_theta(self, self.PS.str2_1{i}, self.PS.str2_1{i}.u, theta);
                
                % Assing values
                self.PS.str2_1{i} = assign_shock_polar(self.PS.str2_1{i}, self.PS.str2{i});
                self.PS.str3_1{i} = assign_shock_polar(self.PS.str3_1{i}, self.PS.strP{i});
                self.PS.str3_2{i} = assign_shock_polar(self.PS.str3_2{i}, self.PS.strP{i});

            else

                try
                    beta = self.PD.beta.value(i);
                catch
                    beta = self.PD.beta.value;
                end

                [self.PS.strR{i}, self.PS.str2{i}] = shock_polar(self, self.PS.strR{i}, u1);
                [~, self.PS.str2_1{i}] = shock_oblique_beta(self, self.PS.strR{i}, u1, beta);
                [~, self.PS.strP{i}] = shock_polar(self, self.PS.str2_1{i}, self.PS.str2_1{i}.u);
                [~, self.PS.str3_1{i}] = shock_oblique_beta(self, self.PS.str2_1{i}, self.PS.str2_1{i}.u, beta);

                % Assing values
                self.PS.str2_1{i} = assign_shock_polar(self.PS.str2_1{i}, self.PS.str2{i});
                self.PS.str3_1{i} = assign_shock_polar(self.PS.str3_1{i}, self.PS.strP{i});
            end

            
        
        case {'SHOCK_POLAR_LIMITRR'}

            try
                u1 = self.PD.u1.value(i);
            catch
                u1 = self.PD.u1.value;
            end
            
            [self.PS.strR{i}, self.PS.str2{i}, self.PS.str2_1{i}, self.PS.strP{i}] = shock_polar_limitRR(self, self.PS.strR{i}, u1);
            % Assing values
            self.PS.str2_1{i} = assign_shock_polar(self.PS.str2_1{i}, self.PS.str2{i});
            
        case {'DET'}

            if i == self.C.l_phi
                [self.PS.strR{i}, self.PS.strP{i}] = det_cj(self, self.PS.strR{i});
            else
                [self.PS.strR{i}, self.PS.strP{i}] = det_cj(self, self.PS.strR{i}, self.PS.strP{i + 1});
            end

        case {'DET_R'}

            if i == self.C.l_phi
                % Calculate post-shock state (2)
                [self.PS.strR{i}, self.PS.str2{i}] = det_cj(self, self.PS.strR{i});
                % Calculate post-shock state (5)
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, self.PS.strR{i}.u, self.PS.str2{i});
            else
                % Calculate post-shock state (2)
                [self.PS.strR{i}, self.PS.str2{i}] = det_cj(self, self.PS.strR{i}, self.PS.str2{i + 1});
                % Calculate post-shock state (5)
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, self.PS.strR{i}.u, self.PS.str2{i}, self.PS.strP{i + 1});
            end

        case {'DET_OVERDRIVEN'}

            try
                drive_factor = self.PD.drive_factor.value(i);
            catch
                drive_factor = self.PD.drive_factor.value;
            end

            if i == self.C.l_phi
                [self.PS.strR{i}, self.PS.strP{i}] = det_overdriven(self, self.PS.strR{i}, drive_factor);
            else
                [self.PS.strR{i}, self.PS.strP{i}] = det_overdriven(self, self.PS.strR{i}, drive_factor, self.PS.strP{i + 1});
            end

        case {'DET_OVERDRIVEN_R'}

            try
                drive_factor = self.PD.drive_factor.value(i);
            catch
                drive_factor = self.PD.drive_factor.value;
            end

            if i == self.C.l_phi
                % Calculate post-shock state (2)
                [self.PS.strR{i}, self.PS.str2{i}] = det_overdriven(self, self.PS.strR{i}, drive_factor);
                % Calculate post-shock state (5)
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, drive_factor, self.PS.str2{i});
            else
                % Calculate post-shock state (2)
                [self.PS.strR{i}, self.PS.str2{i}] = det_overdriven(self, self.PS.strR{i}, drive_factor, self.PS.str2{i + 1});
                % Calculate post-shock state (5)
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, drive_factor, self.PS.str2{i}, self.PS.strP{i + 1});
            end

        case {'DET_UNDERDRIVEN'}

            try
                drive_factor = self.PD.drive_factor.value(i);
            catch
                drive_factor = self.PD.drive_factor.value;
            end

            if i == self.C.l_phi
                [self.PS.strR{i}, self.PS.strP{i}] = det_underdriven(self, self.PS.strR{i}, drive_factor);
            else
                [self.PS.strR{i}, self.PS.strP{i}] = det_underdriven(self, self.PS.strR{i}, drive_factor, self.PS.strP{i + 1});
            end

        case {'DET_UNDERDRIVEN_R'}

            try
                drive_factor = self.PD.drive_factor.value(i);
            catch
                drive_factor = self.PD.drive_factor.value;
            end

            if i == self.C.l_phi
                % Calculate post-shock state (2)
                [self.PS.strR{i}, self.PS.str2{i}] = det_underdriven(self, self.PS.strR{i}, drive_factor);
                % Calculate post-shock state (5)
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, drive_factor, self.PS.str2{i});
            else
                % Calculate post-shock state (2)
                [self.PS.strR{i}, self.PS.str2{i}] = det_underdriven(self, self.PS.strR{i}, drive_factor, self.PS.str2{i + 1});
                % Calculate post-shock state (5)
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, drive_factor, self.PS.str2{i}, self.PS.strP{i + 1});
            end

        case {'DET_OBLIQUE'}

            try
                drive_factor = self.PD.drive_factor.value(i);
            catch
                drive_factor = self.PD.drive_factor.value;
            end

            if ~isempty(self.PD.theta.value)

                try
                    theta = self.PD.theta.value(i);
                catch
                    theta = self.PD.theta.value;
                end

                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = det_oblique_theta(self, self.PS.strR{i}, drive_factor, theta);
            else

                try
                    beta = self.PD.beta.value(i);
                catch
                    beta = self.PD.beta.value;
                end

                [self.PS.strR{i}, self.PS.strP{i}] = det_oblique_beta(self, self.PS.strR{i}, drive_factor, beta);
            end

        case {'DET_POLAR'}

            try
                drive_factor = self.PD.drive_factor.value(i);
            catch
                drive_factor = self.PD.drive_factor.value;
            end

            [self.PS.strR{i}, self.PS.strP{i}] = det_polar(self, self.PS.strR{i}, drive_factor);
        case {'ROCKET'}

            try
                Aratio = self.PD.Aratio.value(i);
            catch
                Aratio = self.PD.Aratio.value;
            end

            if i == self.C.l_phi
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.mix2_c{i}, self.PS.mix3{i}, self.PS.strP{i}] = rocket_performance(self, self.PS.strR{i}, Aratio);
            else
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.mix2_c{i}, self.PS.mix3{i}, self.PS.strP{i}] = rocket_performance(self, self.PS.strR{i}, Aratio, self.PS.str2{i + 1}, self.PS.mix2_c{i + 1}, self.PS.mix3{i + 1}, self.PS.strP{i + 1});
            end

        otherwise

            if length(self.PD.pP.value) > 1
                pP = self.PD.pP.value(i);
            else
                pP = self.PD.pP.value;
            end

            if i == self.C.l_phi
                self.PS.strP{i} = equilibrate(self, self.PS.strR(i), pP);
            else
                self.PS.strP{i} = equilibrate(self, self.PS.strR(i), pP, self.PS.strP(i + 1));
            end

    end

end

function self = check_complete_reaction(self, i)
    % Check if the list of species corresponds to "complete_reaction"
    % If FLAG_COMPLETE is true, establish the list of species based on the
    % given equivalence ratio (phi)
    if self.S.FLAG_COMPLETE
        EquivalenceRatio = self.PS.strR{i}.phi;
        EquivalenceRatio_soot = self.PS.strR{i}.phi_c;

        if EquivalenceRatio < 1
            LS = self.S.LS_lean;
        elseif EquivalenceRatio >= 1 && EquivalenceRatio < EquivalenceRatio_soot
            LS = self.S.LS_rich;
        else
            LS = self.S.LS_soot;
        end

        self.Misc.index_LS_original = find_ind(self.S.LS, LS);
        self = reorganize_index_phase_species(self, LS);
    end

end

function mix = assign_shock_polar(mix, mix_polar)
    % Assign values from the polar curves
    mix.theta_max = mix_polar.theta_max;
    mix.beta_max = mix_polar.beta_max;
    mix.theta_sonic = mix_polar.theta_sonic;
    mix.beta_sonic = mix_polar.beta_sonic;
end
