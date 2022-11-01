classdef eos_PR_classdef
    % Peng-Robinson Equation of State (EoS)
    properties
        Tc          % Critical tempreature [K]
        Pc          % Critical pressure [bar]
        R           % Gas constant [J/K]
        w           % Acentric factor []
        itMax = 30; % Max number of iterations
        tol = 1e-4; % Tolerance Newton-Raphson method
    end

    methods(Static)
        function value = f(v, a, b, T, p, R, alpha)
            value = v^3 + (b - R * T / p) * v^2 - (3*b^2 + 2*R * T * b / p - a * alpha / p) * v + (b^3 + R * T * b^2 - a * b * alpha) / p;
        end

        function value = fprime(v, a, b, T, p, R, alpha)
            value = 3*v^2 + 2*(b - R * T / p) * v - (3*b^2 + 2*R * T * b / p - a * alpha / p);
        end
        
        function value = a(self)
            value = 0.45724 * self.R^2 * self.Tc^2 / self.Pc;
        end

        function value = b(self)
            value = 0.0778 * self.R * self.Tc / self.Pc;
        end

        function value = alpha(self, T)
            Tr = T / self.Tc; % Reduced temperature
            beta = 0.37464 + 1.54226 * self.w - 0.26992 * self.w^2;
            value = (1 + beta * (1 - sqrt(Tr)))^2;
        end
        
        function v_guess = compute_v_molar_guess(T, p, R)
            % Compute molar volume considering ideal EoS
            v_guess = R * T / p;
        end
    end

    methods
        function self = eos_PR(varargin)
            % Unpack inputs
            for i = 1:nargin
                switch lower(varargin{i})
                    case {'tc', 'temperature'}
                        self.Tc = varargin{i+1};
                    case {'pc', 'pressure'}
                        self.Pc = varargin{i+1};
                    case {'r'}
                        self.R = varargin{i+1};
                    case {'w', 'omega', 'acentric factor'}
                        self.w = varargin{i+1};
                    case {'tol', 'tolerance'}
                        self.tol = varargin{i+1};
                    case {'itmax', 'iterations'}
                        self.itMax = varargin{i+1};
                end
            end
        end

        function [a, b, alpha] = get_coefficients(self, T)
            a = self.a(self); % Dimensionless atraction
            b = self.b(self); % Dimensionaless covolume
            alpha = self.alpha(self, T);
        end

        function [f, fprime] = compute_newton_functions(self, v, T, p, R)
            [a, b, alpha] = get_coefficients(self, T);
            f = self.f(v, a, b, T, p, R, alpha);
            fprime = self.fprime(v, a, b, T, p, R, alpha);
        end

        function [x, STOP, it] = compute_v_molar(self, T, p, R)
            x0 = self.compute_v_molar_guess(T, p, R);
            
            it = 0; STOP = 1.0;
            while STOP > self.tol && it < self.itMax
                it = it + 1;
                
                [f0, fprime0] = compute_newton_functions(self, x0, T, p, R);

                x = abs(x0 - f0 / fprime0);
        
                STOP = abs((x - x0) / x);
                x0 = x;
            end
        end
    end
end