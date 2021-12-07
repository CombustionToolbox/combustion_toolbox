function gui_update_from_uitree(obj, selectedNodes)
    % Update GUI with the data that correspond with the selected node from
    % the UITree
    index = regexp(selectedNodes.Text, '[0-9]');
    value = sscanf(selectedNodes.Text(index), '%d');
    results = selectedNodes.Parent.NodeData;
    gui_write_results(obj, results, value, 'FLAG_REACTANTS', true);
end