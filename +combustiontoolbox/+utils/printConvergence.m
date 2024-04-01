function printConvergence(it, itMax, T, STOP, TOL)
    % Print error of the method if the number of iterations is greater than maximum iterations allowed
    %
    % Args:
    %     it (float):    Number of iterations executed in the method
    %     itMax (float): Maximum nNumber of iterations allowed in the method
    %     T (float):     Temperature [K]
    %     STOP (float):  Relative error [-]
    %     tol (float):   Relative tolerance [-]

    if STOP > TOL || it >= itMax
        fprintf('***********************************************************\n')
        fprintf('Root algorithm not converged \n')
        fprintf('   Error       =  %8.2f [%%]  \n', STOP*100)
        fprintf('   Temperature =  %8.2f [K]  \n', T)
        fprintf('   Iterations  =  %8.d [it] \n', it)
        fprintf('***********************************************************\n')
    end

    if T > 2e4
        fprintf('***********************************************************\n')
        fprintf('Validity of the results compromise: T = %d K\n', round(T))
        fprintf('Thermodynamic properties fitted to 20000 K\n');
    end

end