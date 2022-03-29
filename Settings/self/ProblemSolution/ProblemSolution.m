function self = ProblemSolution()
    % Initialize struct with problem solution data
    % 
    % Returns:
    %     self (struct): struct with problem solution data

    % Description
    self.description = "Problem solution";
    % Variables
    self.strR_Fuel = [];     % Properties of the initial fuel mixture             (reactants)
    self.strR_Oxidizer = []; % Properties of the initial oxidizer + inert mixture (reactants)
    self.strR = [];          % Properties of the initial mixture                  (reactants)
    self.strP = [];          % Properties of the final mixture                    (products)
    self.mix1 = [];          % Properties of the mixture at the state 1           (reactants)
    self.mix2 = [];          % Properties of the mixture at the state 2           (products)
    self.mix3 = [];          % Properties of the mixture at the state 3           (products)
    self.mix4 = [];          % Properties of the mixture at the state 4           (products)
end