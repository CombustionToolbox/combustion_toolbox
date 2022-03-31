function varargout = print_stoichiometric_matrix(self, varargin)
    % Print stoichiometric matrix
    %
    % Args:
    %     self (struct): Data of the mixture, conditions, and databases
    %
    % Optional args:
    %     type (str): 'transpose'
    %
    % Optional returns:
    %     A0 (table): Stoichiometric matrix. In case type == 'transpose'
    %                 it returns the transpose of stoichiometric matrix
    
    % Definitions (print)
    label_type = 'Stochiometric matrix:';
    type = ' ';
    % Unpack
    if nargin > 1
        type = varargin{1};
    end
    % Get elements
    elements = self.E.elements';
    % Get species
    species = self.S.LS;
    % Get transpose of the stochiometric matrix
    A0 = self.C.A0.value;
    A0_T = A0';
    % Set table
    T = array2table(A0_T, 'RowNames', elements, 'VariableNames', species);
    % Print table
    if strcmpi(type, 'transpose')
        label_type = 'Transpose stochiometric matrix:' ;
    end
    fprintf('\n%s\n\n', label_type)
    disp(T);
    % Return stochiometric matrix
    varargout = {T};
end
