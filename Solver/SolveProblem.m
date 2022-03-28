function self = SolveProblem(self, ProblemType)
    % Solve the given ProblemType with the conditions and mixture specified in self
    %
    % Args:
    %     self (struct):     Data of the mixture, conditions, and databases
    %     ProblemType (str): Tag of the problem to solve
    %
    % Returns:
    %     self (struct):     Data of the mixtures (initial and final), conditions, databases
    
    try
        % Save Problem Type
        self.PD.ProblemType = ProblemType;
        % Check inputs and set length of the loop
        self = check_inputs(self);
        % Get Flags and length of the loop
        self = get_FLAG_N(self);
        for i = self.C.l_phi:-1:1
            % SET PROBLEM CONDITIONS PER CASE
            self = set_problem_conditions(self, i);
            % COMPUTE PROPERTIES INITIAL MIXTURE
            self = Define_FOI(self, i);
            % CHECK IF LIST OF PRODUCTS CORRESPONDS WITH COMPLETE REACTION
            self = check_complete_reaction(self, i);
            % SOLVE SELECTED PROBLEM
            self = selectProblem(self, i);
            % DISPLAY RESULTS COMMAND WINDOW
            results(self, i);
        end
        % RECOVER RANGE VALUES (ONLY FOR SET OF PROBLEMS)
        if isfield(self.PD, 'range_name')
            self.PD.(self.PD.range_name).value = self.PD.range;
        end
    catch ME
        errorMessage = sprintf('Error in function %s() at line %d.\n\nError Message:\n%s', ...
                               ME.stack(1).name, ME.stack(1).line, ME.message);
        fprintf('%s\n', errorMessage);
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

function self = selectProblem(self, i)
    % Solve selected problem
    switch upper(self.PD.ProblemType)
        case {'SHOCK_I', 'SHOCK_R'}
            try
                u1 = self.PD.u1.value(i);
            catch
                u1 = self.PD.u1.value;
            end
            if strcmp(self.PD.ProblemType,'SHOCK_I')
                if i==self.C.l_phi
                    [self.PS.strR{i}, self.PS.strP{i}] = shock_incident(self, self.PS.strR{i}, u1);
                else
                    [self.PS.strR{i}, self.PS.strP{i}] = shock_incident(self, self.PS.strR{i}, u1, self.PS.strP{i+1});
                end
            else
                if i==self.C.l_phi
                    % Calculate post-shock state (2)
                    [self.PS.strR{i}, self.PS.str2{i}] = shock_incident(self, self.PS.strR{i}, u1);
                    % Calculate post-shock state (5)
                    [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, u1, self.PS.str2{i});
                else
                    % Calculate post-shock state (2)
                    [self.PS.strR{i}, self.PS.str2{i}] = shock_incident(self, self.PS.strR{i}, u1, self.PS.str2{i+1});
                    % Calculate post-shock state (5)
                    [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, u1, self.PS.str2{i}, self.PS.strP{i+1});
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
                [~, self.PS.str2_1{i}, ~] = shock_oblique_theta(self, self.PS.strR{i}, u1, theta);
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
            end
        case {'DET'}
            if i==self.C.l_phi
                [self.PS.strR{i}, self.PS.strP{i}] = cj_detonation(self, self.PS.strR{i});
            else
                [self.PS.strR{i}, self.PS.strP{i}] = cj_detonation(self, self.PS.strR{i}, self.PS.strP{i+1});
            end
        case {'DET_R'}
            if i==self.C.l_phi
                % Calculate post-shock state (2)
                [self.PS.strR{i}, self.PS.str2{i}] = cj_detonation(self, self.PS.strR{i});
                % Calculate post-shock state (5)
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, self.PS.strR{i}.u, self.PS.str2{i});
            else
                % Calculate post-shock state (2)
                [self.PS.strR{i}, self.PS.str2{i}] = cj_detonation(self, self.PS.strR{i}, self.PS.str2{i+1});
                % Calculate post-shock state (5)
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, self.PS.strR{i}.u, self.PS.str2{i}, self.PS.strP{i+1});
            end
        case 'DET_OVERDRIVEN'
            try
                overdriven = self.PD.overdriven.value(i);
            catch
                overdriven = self.PD.overdriven.value;
            end
            if i==self.C.l_phi
                [self.PS.strR{i}, self.PS.strP{i}] = overdriven_detonation(self, self.PS.strR{i}, overdriven);
            else
                [self.PS.strR{i}, self.PS.strP{i}] = overdriven_detonation(self, self.PS.strR{i}, overdriven, self.PS.strP{i+1});
            end
        case {'DET_OBLIQUE'}
            try
                overdriven = self.PD.overdriven.value(i);
            catch
                overdriven = self.PD.overdriven.value;
            end
            if ~isempty(self.PD.theta.value)
                try
                    theta = self.PD.theta.value(i);
                catch
                    theta = self.PD.theta.value;
                end
                    [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = det_oblique_theta(self, self.PS.strR{i}, overdriven, theta);
            else
                try
                    beta = self.PD.beta.value(i);
                catch
                    beta = self.PD.beta.value;
                end
                [self.PS.strR{i}, self.PS.strP{i}] = det_oblique_beta(self, self.PS.strR{i}, overdriven, beta);
            end
        case 'DET_OVERDRIVEN_R'
            try
                overdriven = self.PD.overdriven.value(i);
            catch
                overdriven = self.PD.overdriven.value;
            end
            if i==self.C.l_phi
                % Calculate post-shock state (2)
                [self.PS.strR{i}, self.PS.str2{i}] = overdriven_detonation(self, self.PS.strR{i}, overdriven);
                % Calculate post-shock state (5)
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, overdriven, self.PS.str2{i});
            else
                % Calculate post-shock state (2)
                [self.PS.strR{i}, self.PS.str2{i}] = overdriven_detonation(self, self.PS.strR{i}, overdriven, self.PS.str2{i+1});
                % Calculate post-shock state (5)
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = shock_reflected(self, self.PS.strR{i}, overdriven, self.PS.str2{i}, self.PS.strP{i+1});
            end
        case 'ROCKET'
            if i==self.C.l_phi
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = rocket_performance(self, self.PS.strR{i});
            else
                [self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i}] = rocket_performance(self, self.PS.strR{i}, self.PS.str2{i+1}, self.PS.strP{i+1});
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
    mix.theta_max = mix_polar.theta_max;
    mix.beta_max =  mix_polar.beta_max;
    mix.theta_sonic = mix_polar.theta_sonic;
    mix.beta_sonic = mix_polar.beta_sonic;
end
