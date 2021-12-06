function gui_addchildren(app, app_parent, results)
    %%% New Children Node: Results -> ProblemType
    ProblemType = results.ProblemType;
    % Category 1
    j = 0;
    if ~isempty(app_parent.Children)
        i = 1;
        while i<= length(app_parent.Children)
            if strcmpi(app_parent.Children(i).Text, ProblemType)
                category1 = app_parent.Children(i);
                j = 1;
                break;
            end
            i = i + 1;
        end
    end
    if j == 0
        category1 = uitreenode(app_parent, 'Text', ProblemType, 'NodeData', []);
    end
    expand(app_parent,'all');
    %%% New Children Node: ProblemType -> Reactants
    for k=results.length:-1:1
        % Category 2
        j = 0;
        text_reactants = app.Reactants.Items{sscanf(results(k).numberReactants, '%d')};
        if ~isempty(category1.Children)
            i = 1;
            while i<=length(category1.Children)
                if strcmpi(category1.Children(i).Text, text_reactants)
                    category2 = category1.Children(i);
                    j = 1;
                    break;
                end
                i = i + 1;
            end
        end
        if j == 0
            category2 = uitreenode(category1, 'Text', text_reactants, 'NodeData', []);
        end
        expand(category1,'all');
        %%% New Children Node: Reactants -> CompleteOrIncomplete
        % Category 3
        j = 0;
        if ~isempty(category2.Children)
            i = 1;
            while i<=length(category2.Children)
                if strcmpi(category2.Children(i).Text, results.reaction)
                    category3 = category2.Children(i);
                    j = 1;
                    break;
                end
                i = i + 1;
            end
        end
        if j == 0
            category3 = uitreenode(category2, 'Text', results.reaction, 'NodeData', [], 'UserData', 1);
        end
        expand(category2,'all');
        %%% New Children Node: CompleteOrIncomplete -> Phi
        % Category 4
%                 j = 0;
        %             if any(strcmpi(results.app.NR,{'1','2'})) % Cases without fuel
        text = 'Case %d'; text_value = category3.UserData;
        %             end
%                 if ~isempty(category3.Children)
%                     i = 1;
%                     while i<=length(category3.Children)
%                         if strcmpi(category3.Children(i).Text,sprintf(text,text_value))
% %                             category4 = category3.Children(i);
%                             j = 1;
%                             break;
%                         end
%                         i = i + 1;
%                     end
%                 end
%                 if j == 0
            uitreenode(category3, 'Text', sprintf(text,text_value), 'NodeData', results);
            category3.UserData = category3.UserData + 1;
%                 end
        expand(category3,'all');
    end
end