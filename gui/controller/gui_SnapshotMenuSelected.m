function gui_SnapshotMenuSelected(UIFigure)
    filter = {'*.pdf';'*.jpg';'*.png';'*.tif'};
    [filename, filepath] = uiputfile(filter);
    if ischar(filename)
        exportapp(UIFigure, [filepath filename]);
    end
end