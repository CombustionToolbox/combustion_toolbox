function g0 = get_g0(app, species, T, DB)
    if DB.(species).ctTInt > 0
        R0 = 8.3144598;
        a = DB.(species).a;
        b = DB.(species).b;
        ctTInt       = DB.(species).ctTInt;
        tRange       = DB.(species).tRange;
        tExponents   = DB.(species).tExponents;
        Element_matrix = set_element_matrix(DB.(species).txFormula,app.E.Elements);
        
        if (T < tRange{1}(1)) || (T > tRange{ctTInt}(2))
            cP = species_cP(species, T, DB);
            h0 = species_DhT(species, T, DB);
            s0 = species_s0(species, T, DB);
            g0 = species_g0(species, T, DB);
            return
        else
            for i = 1:ctTInt
                if (T >= tRange{i}(1)) && (T <= tRange{i}(2))
                    tInterval = i;
                end
            end
            cP = R0 *      sum(a{tInterval} .* T.^tExponents{tInterval});
            h0 = R0 * T * (sum(a{tInterval} .* T.^tExponents{tInterval} .* [-1   log(T) 1      1/2 1/3 1/4 1/5 0]) + b{tInterval}(1)/T);
            s0 = R0 *     (sum(a{tInterval} .* T.^tExponents{tInterval} .* [-1/2 -1     log(T) 1   1/2 1/3 1/4 0]) + b{tInterval}(2)  );
            gP = h0 - T.*s0;
            for i = 1:size(Element_matrix,2)
            nu_i = Element_matrix(2,i);
            [iRE_i, REname_i] = isRefElm(Reference_form_of_elements_with_T_intervals,Elements{Element_matrix(1,i)},T);
%             [~, REname_i] = isRefElm(Reference_form_of_elements_with_T_intervals,upper(Elements{Element_matrix(1,i)}),T);
%             [~, ~, ~, ~, ~, H0_i, ~, ~, S0_i, ~] = SpeciesThermProp(strMaster,REname_i,T,'molar',0);
            [txFormula_i, mm_i, Cp0_i, Cv0_i, Hf0_i, H0_i, Ef0_i, E0_i, S0_i, DfG0_i] = SpeciesThermProp(strMaster,REname_i,T,'molar',0);
            GR(i) = nu_i*(H0_i - T.*S0_i);
            if any(Element_matrix(1,i)==[1, 7, 8, 9, 17, 35]), GR(i) = GR(i)/2; end
            end
            GR = sum(GR);
            DfG0 = GP-GR;
            
            h0 = h0*1e-3;
            g0 = g0*1e-3;
        end
    end
end
