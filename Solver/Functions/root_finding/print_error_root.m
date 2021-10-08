function print_error_root(it, itMax, TP, ERR)
    if it == itMax
        fprintf('***********************************************************\n')
        fprintf('Root algorithm not converged \n')
        fprintf('   Error       =  %8.2f [%%]  \n', ERR*100)
        fprintf('   Temperature =  %8.2f [K]  \n', TP)
        fprintf('   Iterations  =  %8.d [it] \n', it)
        fprintf('***********************************************************\n')
    end
end