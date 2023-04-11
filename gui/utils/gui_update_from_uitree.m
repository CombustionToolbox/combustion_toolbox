function gui_update_from_uitree(app, selectedNodes)
    % Update GUI with the data that correspond with the selected node from
    % the UITree

    % Check if the selected node has children (i.e., does not contain results
    % from a previous calculation)
    if ~isempty(selectedNodes.Children)
        return
    end

    % Get the results from the selected node
    results = selectedNodes.NodeData;

    % Update the GUI with the results
    gui_write_results(app, results, 1, 'FLAG_REACTANTS', true);
end