function gui_SnapshotMenuSelected(UIFigure)
    % Routine to exports the current figure to a file. This
    % function is called when the Snapshot menu is selected.
    % 
    % Notes:
    %      * The file type is determined by the file extension
    %      * The file name is determined by the user
    %      * The file path is determined by the user

    % Definitions
    filter = {'*.pdf';'*.jpg';'*.png';'*.tif'};

    % Get the file name and path
    [filename, filepath] = uiputfile(filter);

    % Export the figure if the user did not cancel
    if ischar(filename)
        exportapp(UIFigure, [filepath filename]);
    end

end