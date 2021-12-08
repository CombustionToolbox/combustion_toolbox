function gui_compute_mach_or_velocity(obj, inputname)
    % Function that computes the pre-shock Mach number of the mixture from 
    % a given pre-shock velocity or viceversa.
    if any(contains(obj.ProblemType.Value, 'SHOCK', 'IgnoreCase', true))
        temp_app = gui_create_temp_app(obj, [], false);
        mix1 = temp_app.PS.strR{1};

        if strcmpi(inputname, 'Mach')
            M1 = gui_get_prop(obj, 'M1', obj.PR4.Value);
            u1 = M1 * soundspeed(mix1);
            obj.PR3.Value = sprintf('%.6g', u1);
        else
            u1 = gui_get_prop(obj, 'u1', obj.PR3.Value);
            M1 = u1 / soundspeed(mix1);
            obj.PR4.Value = sprintf('%.6g', M1);
        end
    end
end