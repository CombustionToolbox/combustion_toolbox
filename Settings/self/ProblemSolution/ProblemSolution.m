function self = ProblemSolution()
    % Initialize struct with problem solution data
    % 
    % Returns:
    %     self (struct): struct with problem solution data

    % Description
    self.description = "Problem solution";
    % Variables
    self.strR_Fuel = [];               % Properties of the initial fuel mixture             (reactants)
    self.strR_Oxidizer = [];           % Properties of the initial oxidizer + inert mixture (reactants)
    self.strR = [];                    % Properties of the initial mixture                  (reactants)
    self.strP = [];                    % Properties of the final mixture                    (products)
end