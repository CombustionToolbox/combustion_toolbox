function print_error_root(it, itMax, TP, ERR)
    if it > itMax
        print('**********************************\n')
        print('** Root algorithm not converged **\n')
        print('** Temp  =  %4.2f               **\n', TP)
        print('** Error =  %4.2f%%             **\n', abs(ERR)*100)
        print('** It    =  %4.d                **\n', it)
        print('**********************************\n')
    end
end