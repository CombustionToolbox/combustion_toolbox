function print_mixture(self)
    fprintf('\n\nMixture:\n\n');
    fprintf('temperature       |\t   %8.3f  %s\n', self.T,     'K');
    fprintf('pressure          |\t   %8.3f  %s\n', self.p,     'bar');
    fprintf('density           |\t   %8.3f  %s\n', self.rho,   'kg/m3');
    fprintf('mean mol. weight  |\t   %8.3f  %s\n', self.W,     'amu');
    fprintf('adibatic index    |\t   %8.3f  %s\n', self.gamma, '-');
    fprintf('\n\n');
    fprintf('                  |\t   1 kg             1 mol   \n');
    fprintf('                  |\t-----------      -----------\n');
    fprintf('enthalpy          |\t   %8.3f         %8.3f  %s\n', self.h/self.mi, self.h, 'kJ');
    fprintf('internal energy   |\t   %8.3f         %8.3f  %s\n', self.e/self.mi, self.e, 'kJ');
    fprintf('entropy           |\t   %8.3f         %8.3f  %s\n', self.S, self.S*self.mi, 'kJ/K');
    fprintf('heat capacity c_p |\t   %8.3f         %8.3f  %s\n', self.cP/self.mi*1e-3, self.cP*1e-3, 'kJ/K');
    fprintf('heat capacity c_v |\t   %8.3f         %8.3f  %s\n', self.cV/self.mi*1e-3, self.cV*1e-3, 'kJ/K');  
end

function print_fractions(self)
    %%%% SORT SPECIES COMPOSITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [self.Xi(:),ind_sort] = sort(self.Xi(:), 'descend');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    j = self.Xi>mintol_display;
    Xminor = sum(self.Xi(~j));
    for i=1:length(j)
        if j(i)
            fprintf('%-12s \t%11.4e\n',NameSpecies{ind_sort(i)},self.Xi(i));
%             fprintf('%-10s \t%11.4e\n',NameSpecies{i},strR.Xi(i))
        end
    end
    fprintf('MINORS[+%d]   %12.4e\n\n',length(self.Xi)-sum(j),Xminor);
    fprintf('TOTAL  \t\t %14.4e\n',sum(self.Xi));
end