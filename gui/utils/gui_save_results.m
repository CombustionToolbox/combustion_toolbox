function gui_save_results(app, format)
    % Save results as the given format (.xls or .mat)
    %
    % Args:
    %     app (obj): class with all the data of the app
    %     format (char): format extension/name
    
    % Check if there is data available to export
    if isempty(app.Node_Results.Children)
        app.Console_text.Value = 'Nothing to export.';
        return
    end
    
    % Read format
    switch lower(format)
        case {'xls', '.xls', 'excel'}
            FLAG_EXCEL = true;
            format = '*.xls';
        case {'.mat', 'mat', 'matlab'}
            FLAG_EXCEL = false;
            format = '*.mat';
    end
    
    % Organize data and export to a spreadsheet if FLAG_EXCEL is true
    [data_mix1, data_mix2, filename, FLAG_CANCEL] = gui_save_format_results(app, format, FLAG_EXCEL);
    
    if FLAG_CANCEL
        return
    end

    if FLAG_EXCEL
        return
    end
    
    % Save data in a .mat file
    gui_save_mat(filename, data_mix1, data_mix2)
end

% SUB-PASS FUNCTIONS
 function [data_mix1, data_mix2, filename, FLAG_CANCEL] = gui_save_format_results(app, format, FLAG_EXCEL)
    % Save results as the given format (.xls or .mat)
    
    % Initialization
    data_mix1 = [];
    data_mix2 = [];
    filename = [];
    FLAG_CANCEL = false;

    % Read data
    results = app.Tree.Children;
    sumcases = 0;
    if isempty(results.Children)
        return
    end
    
    % Get child1
    child1 = results.Children;

    % Ask for filename
    [file, name, path] = uiputfile(format, 'File Selection', app.Misc.export_results.filename);
    
    % Check if the users has cancelled the export
    if isfloat(file)
        FLAG_CANCEL = true;
        app.Console_text.Value = 'Export cancelled.';
        return
    end

    filename = [name, file];

    disp(['User selected ', fullfile(path, file), ' and then clicked Save.'])

    % Gather datasets and export them as a spreadsheet
    for i = length(child1):-1:1
        % Get child2
        child2 = child1(i).Children;

        for j = length(child2):-1:1
            % Get child3
            child3 = child2(j).Children;

            for k = length(child3):-1:1
                % Get child4
                child4 = child3(k).Children;

                % Get list of species
                LS = child4(1).NodeData.LS;
                
                % Get mixtures and equivalence ratio
                for l = length(child4):-1:1
                    mix1{l} = child4(l).NodeData.mix1;
                    mix2{l} = child4(l).NodeData.mix2;
                    phi(l)  = child4(l).NodeData.mix1.phi;
                end
                
                % Update sheet number
                sumcases = sumcases + 1;

                % Get data in cell format
                data_mix1{sumcases} = get_excel_cell(mix1, LS, phi);
                data_mix2{sumcases} = get_excel_cell(mix2, LS, phi);

                % Assign sheet id
                sheet_name = sprintf('%s-%d_%d', child1(i).Text(14:end), j, k);
                
                if ~FLAG_EXCEL
                    continue
                end
                
                % Save data in a spreadsheet
                gui_save_excel(filename, data_mix1, data_mix2, sumcases, sheet_name)

            end

        end

    end

end

function gui_save_mat(filename, data_mix1, data_mix2)
    % Save results as a mat file
    save(filename, 'data_mix1', 'data_mix2');
end

function gui_save_excel(filename, data_mix1, data_mix2, sumcases, sheet_name)
    % Save results as an excel file
    writecell(data_mix1{sumcases}, filename, 'Sheet', strcat('R-', sheet_name));
    writecell(data_mix2{sumcases}, filename, 'Sheet', strcat('P-', sheet_name));
end