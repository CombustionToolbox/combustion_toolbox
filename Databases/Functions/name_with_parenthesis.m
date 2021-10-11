function name_with = name_with_parenthesis(name_wo)

    name_i = name_wo;
    j = find(name_i(:)=='b');
    if length(j)==2
        name_i(j(1)) = '(';
        name_i(j(2)) = ')';
    end
    if (length(j)==3)
        if (j(3)-j(1)==2)||((j(2)-j(1)==1)&&(j(1)>2))
           name_i(j(1)) = '(';
           name_i(j(3)) = ')';
        elseif (j(2)-j(1)==2)&&(j(1)>2)       
           name_i(j(1)) = '(';
           name_i(j(2)) = ')';
        else
           name_i(j(2)) = '(';
           name_i(j(3)) = ')';
        end
    end
    if (length(j)==4)
        if j(4)-j(2)==2
            name_i(j(2)) = '(';
            name_i(j(4)) = ')';
        else
            name_i(j(1)) = '(';
            name_i(j(2)) = ')';
            name_i(j(3)) = '(';
            name_i(j(4)) = ')';
        end
    end
    name_with = name_i;