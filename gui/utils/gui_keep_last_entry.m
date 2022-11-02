function gui_keep_last_entry(obj, entry_obj)
    % Function that keeps only the last entry, i.e., if the order is
    % entry1 and entry2, entry1 will be set to empty value

    if strcmpi(obj.ProblemType.Value, 'ROCKET')
        entry_obj.Value = '';
    end
    
end