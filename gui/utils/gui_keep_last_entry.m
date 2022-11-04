function gui_keep_last_entry(app, entry_app)
    % Function that keeps only the last entry, i.e., if the order is
    % entry1 and entry2, entry1 will be set to empty value

    if strcmpi(app.ProblemType.Value, 'ROCKET')
        entry_app.Value = '';
    end
    
end