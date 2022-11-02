function results(self, i)
    % Display results in the command window

    if self.Misc.FLAG_RESULTS

        if ~isfield(self.PS, 'str2')
            print_mixture(self, self.PS.strR{i}, self.PS.strP{i});
        elseif isfield(self.PS, 'mix2_c') && isempty(self.PS.strP{i}) && self.PD.FLAG_IAC
            print_mixture(self, self.PS.strR{i}, self.PS.mix2_c{i}, self.PS.mix3{i});
        elseif isfield(self.PS, 'mix2_c') && ~isempty(self.PS.strP{i}) && self.PD.FLAG_IAC
            print_mixture(self, self.PS.strR{i}, self.PS.mix2_c{i}, self.PS.mix3{i}, self.PS.strP{i});
        elseif isfield(self.PS, 'mix2_c') && ~isempty(self.PS.strP{i}) && ~self.PD.FLAG_IAC
            print_mixture(self, self.PS.strR{i}, self.PS.str2{i}, self.PS.mix2_c{i}, self.PS.mix3{i}, self.PS.strP{i});
        elseif isfield(self.PS, 'mix2_c') && isempty(self.PS.strP{i}) && ~self.PD.FLAG_IAC
            print_mixture(self, self.PS.strR{i}, self.PS.str2{i}, self.PS.mix2_c{i}, self.PS.mix3{i});
        elseif ~isfield(self.PS, 'str3_2')
            print_mixture(self, self.PS.strR{i}, self.PS.str2{i}, self.PS.strP{i});
        elseif ~isfield(self.PS, 'str3')
            print_mixture(self, self.PS.strR{i}, self.PS.str2_1{i}, self.PS.str3_1{i}, self.PS.str3_2{i});
        else
            print_mixture(self, self.PS.strR{i}, self.PS.str2{i}, self.PS.str3{i}, self.PS.str3_2{i});
        end

    end

end
