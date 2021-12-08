function gui_compute_mach_or_velocity(obj, inputname)
    % Function that computes the pre-shock Mach number of the mixture from 
    % a given pre-shock velocity or viceversa.
    if any(contains(obj.ProblemType.Value, 'SHOCK', 'IgnoreCase', true))
        if ~isempty(obj.UITable_R.Data)
            temp_app = gui_create_temp_app(obj, [], false);
            mix1 = temp_app.PS.strR{1};
    
            if strcmpi(inputname, 'Mach')
                [M1, FLAG] = gui_get_prop(obj, 'M1', obj.PR4.Value);
                u1 = M1 * soundspeed(mix1);
                obj.PR3.Value = compute_vector_or_scalar(u1, FLAG);
            else
                [u1, FLAG] = gui_get_prop(obj, 'u1', obj.PR3.Value);
                M1 = u1 / soundspeed(mix1);
                obj.PR4.Value = compute_vector_or_scalar(M1, FLAG);
            end
        end
    end
end

% SUB-PASS FUNCTIONS
function value_char = compute_vector_or_scalar(value, FLAG)
    if FLAG
        value_char = sprintf('[%.5g:%.5g:%.5g]', value(1), value(2) - value(1), value(end));
    else
        value_char = sprintf('%.5g', value);
    end
end