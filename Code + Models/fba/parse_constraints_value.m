
function [constraints_value_type,constraints_value_value,constraints_value_reaction,constraints_value_count] = parse_constraints_value(input_file)
  
    % load data
    [~,~,raw] = xlsread(input_file,'Constraints - Value');

    % value constraints
    constraints_value_type = {};
    constraints_value_value = {};
    constraints_value_reaction = {};
    constraints_value_count = 0;

    for row = 3:size(raw,1)
        added = false;
        
        % constraint type - equal
        value = raw(row,1);
        value = value{1};
        if ~isnan(value)
            if isfloat(value)           
                if value >= 0
                    added = true;
                    constraints_value_type{end+1} = 'equal';
                    constraints_value_value{end+1} = value;
                else
                    error('ERROR - Constraints - Value - Constraint value must be >= 0')
                end
            else
                error('ERROR - Constraints - Value - Constraint value must be a number')
            end
        end

        % constraint type - less than
        value = raw(row,2);
        value = value{1};
        if ~isnan(value)
            if ~added;
                if isfloat(value)           
                    if value >= 0
                        added = true;
                        constraints_value_type{end+1} = 'less';
                        constraints_value_value{end+1} = value;
                    else
                        error('ERROR - Constraints - Value - Constraint value must be >= 0')
                    end
                else
                    error('ERROR - Constraints - Value - Constraint value must be a number')
                end
            else
                error('ERROR - Constraints - Value - Can only choose one constraint type at a time')
            end
        end

        % constraint type - greater than
        value = raw(row,3);
        value = value{1};
        if ~isnan(value)
            if ~added;
                if isfloat(value)           
                    if value >= 0
                        added = true;
                        constraints_value_type{end+1} = 'greater';
                        constraints_value_value{end+1} = value;
                    else
                        error('ERROR - Constraints - Value - Constraint value must be >= 0')
                    end
                else
                    error('ERROR - Constraints - Value - Constraint value must be a number')
                end
            else
                error('ERROR - Constraints - Value - Can only choose one constraint type at a time')
            end
        end

        % constraint type - range
        value1 = raw(row,4);
        value1 = value1{1};
        value2 = raw(row,5);
        value2 = value2{1};
        if ~isnan(value1) && ~isnan(value2)
            if ~added;
                if isfloat(value1) && isfloat(value2)       
                    if (value1 >= 0) && (value2 >= 0)
                        if value1 <= value2
                            added = true;
                            constraints_value_type{end+1} = 'range';
                            constraints_value_value{end+1} = [value1,value2];
                        else
                            error('ERROR - Constraints - Value - Minimum value must be <= maximum value')
                        end
                    else
                        error('ERROR - Constraints - Value - Constraint values must be >= 0')
                    end
                else
                    error('ERROR - Constraints - Value - Constraint values must be a number')
                end
            else
                error('ERROR - Constraints - Value - Can only choose one constraint type at a time')
            end
        elseif ~isnan(value1) || ~isnan(value2)
            error('ERROR - Constraints - Value - If choosing "Range" constraint, must enter both minimum and maximum values');
        end

        % reaction formula
        if added
            value = raw(row,6);
            value = value{1};
            if ~isnan(value)
                if ischar(value)
                    constraints_value_reaction{end+1} = value;
                else
                    error('ERROR - Constraints - Value - Reaction formula must be valid string')
                end
            else
                error('ERROR - Constraints - Value - Must give reaction formula');
            end
        end
        
        % constraints counter
        if added
            constraints_value_count = constraints_value_count + 1;
        end
    end
end