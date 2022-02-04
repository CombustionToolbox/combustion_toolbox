function gui_save_results(app, format)
    % Save results as the given format (.xls or .mat)

    switch lower(format)
        case {'xls', '.xls', 'excel'}
            FLAG_EXCEL = true;
            format = '*.xls';
        case {'.mat', 'mat', 'matlab'}
            FLAG_EXCEL = false;
            format = '*.mat';
    end

    [data_mix1, data_mix2, filename] = gui_save_format_results(app, format, FLAG_EXCEL);
    
    if FLAG_EXCEL
        return
    end

    gui_save_mat(filename, data_mix1, data_mix2)
end

% SUB-PASS FUNCTIONS
 function [data_mix1, data_mix2, filename] = gui_save_format_results(app, format, FLAG_EXCEL)
    results = app.Tree.Children;
    sumcases = 0;
    if ~isempty(results.Children)
        child1 = results.Children;
        maxi = numel(child1);
        [file, name, path] = uiputfile(format, 'File Selection', '');
        filename = strcat(name, file);
        if isequal(file,0) || isequal(path,0)
            disp('User clicked Cancel.')
        else
            disp(['User selected ', fullfile(path,file), ' and then clicked Save.'])
            for i = maxi:-1:1
                child2 = child1(i).Children;
                for j=numel(child2):-1:1
                    child3 = child2(j).Children;
                    for k=numel(child3):-1:1
                        child4 = child3(k).Children;
                        LS = child4(1).NodeData.LS;
                        for l=numel(child4):-1:1
                            mix1{l} = child4(l).NodeData.mix1;
                            mix2{l} = child4(l).NodeData.mix2;
                            phi(l)  = child4(l).NodeData.mix1.phi;
                        end
                        sumcases = sumcases + 1;
                        data_mix1{sumcases} = FormattedOutput_test([], phi, mix1, LS);
                        data_mix2{sumcases} = FormattedOutput_test([], phi, mix2, LS);
                        sheet_name = strcat(child1(i).Text(14:end), '-', sprintf('%d', j), '_', sprintf('%d', k));
                        if FLAG_EXCEL
                            gui_save_excel(filename, data_mix1, data_mix2, sumcases, sheet_name)
                        end
                    end
                end
            end
        end
    end
end

function gui_save_mat(filename, data_mix1, data_mix_2)
    % Save results as a mat file
    save(filename, 'data_mix1', 'data_mix_2');
end

function gui_save_excel(filename, data_mix1, data_mix2, sumcases, sheet_name)
    % Save results as an excel file
    writecell(data_mix1{sumcases}, filename, 'Sheet', strcat('R-', sheet_name));
    writecell(data_mix2{sumcases}, filename, 'Sheet', strcat('P-', sheet_name));
end