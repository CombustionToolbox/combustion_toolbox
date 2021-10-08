function print_error_root(it, itMax, TP, ERR)
    if it == itMax
        fprintf('***********************************************************\n')
        fprintf('Root algorithm not converged \n')
        fprintf('   Temp  =  %8.2f            \n', TP)
        fprintf('   Error =  %8.2f            \n', abs(ERR)*100)
        fprintf('   It    =  %8.d             \n', it)
        fprintf('***********************************************************\n')
    end
end