classdef eos_VW
    % Van der Waals Equation of State (EoS)
    properties
        Tc
        Pc
        R
        itMax = 30;
        tol = 1e-4;
    end

    methods(Static)
        function value = f(v, a, b, T, p, R)
            value = v^3 - (b + R * T / p) * v^2 + (a / p) * v - a * b / p;
        end

        function value = fprime(v, a, b, T, p, R)
            value = 3*v^2 - 2*(b + R * T / p) * v + (a / p);
        end
        
        function value = a(self)
            value = 27/64 * self.R^2 * self.Tc^2 / self.Pc;
        end

        function value = b(self)
            value = self.R * self.Tc / (8 * self.Pc);
        end

        function v_guess = compute_v_molar_guess(T, p, R)
            v_guess = R * T / p;
        end
    end

    methods
        function self = eos_VW(varargin)
            % Unpack inputs
            for i = 1:nargin
                switch lower(varargin{i})
                    case {'tc', 'temperature'}
                        self.Tc = varargin{i+1};
                    case {'pc', 'pressure'}
                        self.Pc = varargin{i+1};
                    case {'r'}
                        self.R = varargin{i+1};
                    case {'tol', 'tolerance'}
                        self.tol = varargin{i+1};
                    case {'itmax', 'iterations'}
                        self.itMax = varargin{i+1};
                end
            end
        end

        function [a, b] = get_coefficients(self)
            a = self.a(self);
            b = self.b(self);
        end

        function [f, fprime] = compute_newton_functions(self, v, T, p, R)
            self.R = R;
            [a, b] = get_coefficients(self);
            f = self.f(v, a, b, T, p, R);
            fprime = self.fprime(v, a, b, T, p, R);
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