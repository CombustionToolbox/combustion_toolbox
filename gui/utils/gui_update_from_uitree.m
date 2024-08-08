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

    % Update GUI layout if required
    switch upper(results.mix2.problemType)
        case 'ROCKET_IAC'
            app.ProblemType.Value = 'ROCKET';
            app.FLAG_IAC.Value = true;
        case 'ROCKET_FAC'
            app.ProblemType.Value = 'ROCKET';
            app.FLAG_IAC.Value = false;
        otherwise
            app.ProblemType.Value = results.mix2.problemType;
    end

    gui_ProblemTypeValueChanged(app);

    % Update the GUI with the results
    gui_write_results(app, results, 1, 'FLAG_REACTANTS', true);
end