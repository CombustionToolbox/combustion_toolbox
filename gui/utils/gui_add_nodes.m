function gui_add_nodes(parent_node, results)
    % Function that generate nodes in a UITree and save data on them
    %
    % Args:
    %     parent_node (object): Parent node of the UITree
    %     results (struct): struct with the data to save in the nodes

    % * Category 1: new children node: Results -> ProblemType
    match_text = strcat('Problem Type: ', results(1).ProblemType);
    category1 = set_node(parent_node, match_text);

    % * Category 2: new children node: ProblemType -> Reactants
    match_text = strcat('Reactants: ', results(1).Reactants);
    category2 = set_node(category1, match_text);

    % * Category 3: new children node: Reactants -> ListProducts
    match_text = strcat('List Products: ', results(1).Products);
    category3 = set_node(category2, match_text);

    % * Category 4: new children node: ListProducts -> cases
    for i = 1:length(results)
        match_text = sprintf('Case %d', category3.UserData); 
        set_node(category3, match_text, 'NodeData', results(i));
        category3.UserData = category3.UserData + 1;
    end

end

% SUB-PASS FUNCTIONS
function category_children = set_node(category_parent, match_text, varargin)
    % Default values
    NodeData = [];
    FLAG_CHECK = true;
    % Check varargin inputs
    for i = 1:nargin - 3
        switch lower(varargin{i})
            case 'nodedata'
                NodeData = varargin{i + 1};
                FLAG_CHECK = false;
        end

    end

    % Set node
    j = 0;
    if ~isempty(category_parent.Children) && FLAG_CHECK
        i = 1;
        while i <= length(category_parent.Children)
            if strcmpi(category_parent.Children(i).Text, match_text)
                category_children = category_parent.Children(i);
                j = 1;
                break;
            end

            i = i + 1;
        end

    end

    if j == 0
        category_children = uitreenode(category_parent, 'Text', match_text, 'NodeData', NodeData, 'UserData', 1);
    end

    expand(category_parent,'all');
end