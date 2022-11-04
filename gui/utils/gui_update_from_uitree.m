function gui_update_from_uitree(app, selectedNodes)
    % Update GUI with the data that correspond with the selected node from
    % the UITree
    if isempty(selectedNodes.Children)
        results = selectedNodes.NodeData;
        gui_write_results(app, results, 1, 'FLAG_REACTANTS', true);
    end
end