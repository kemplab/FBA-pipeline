
function [constraints_fraction_type,constraints_fraction_value,constraints_fraction_reaction,constraints_fraction_count] = parse_constraints_fraction(input_file)
  
    % load data
    [~,~,raw] = xlsread(input_file,'Constraints - Fraction');

    % fraction of maximum constraints
    constraints_fraction_type = {};
    constraints_fraction_value = {};
    constraints_fraction_reaction = {};
    constraints_fraction_count = 0;

    for row = 3:size(raw,1)
        added = false;

        % constraint type - equal
        value = raw(row,1);
        value = value{1};
        if ~isnan(value)
            if isfloat(value)           
                if (value >= 0) && (value <= 1)
                    added = true;
                    constraints_fraction_type{end+1} = 'equal';
                    constraints_fraction_value{end+1} = value;
                else
                    error('ERROR - Constraints - Fraction - Constraint fraction must be >= 0 and <= 1')
                end
            else
                error('ERROR - Constraints - Fraction - Constraint fraction must be a number')
            end
        end

        % constraint type - less than
        value = raw(row,2);
        value = value{1};
        if ~isnan(value)
            if ~added
                if isfloat(value)           
                    if (value >= 0) && (value <= 1)
                        added = true;
                        constraints_fraction_type{end+1} = 'less';
                        constraints_fraction_value{end+1} = value;
                    else
                        error('ERROR - Constraints - Fraction - Constraint fraction must be >= 0 and <= 1')
                    end
                else
                    error('ERROR - Constraints - Fraction - Constraint fraction must be a number')
                end
            else
                error('ERROR - Constraints - Fraction - Can only choose one constraint type at a time')
            end
        end

        % constraint type - greater than
        value = raw(row,3);
        value = value{1};
        if ~isnan(value)
            if ~added
                if isfloat(value)           
                    if (value >= 0) && (value <= 1)
                        added = true;
                        constraints_fraction_type{end+1} = 'greater';
                        constraints_fraction_value{end+1} = value;
                    else
                        error('ERROR - Constraints - Fraction - Constraint fraction must be >= 0 and <= 1')
                    end
                else
                    error('ERROR - Constraints - Fraction - Constraint fraction must be a number')
                end
            else
                error('ERROR - Constraints - Fraction - Can only choose one constraint type at a time')
            end
        end

        % constraint type - range
        value1 = raw(row,4);
        value1 = value1{1};
        value2 = raw(row,5);
        value2 = value2{1};
        if ~isnan(value1) && ~isnan(value2)
            if ~added
                if isfloat(value1) && isfloat(value2)       
                    if (value1 >= 0) && (value1 <= 1) && (value2 >= 0) && (value2 <= 1)
                        if value1 <= value2
                            added = true;
                            constraints_fraction_type{end+1} = 'range';
                            constraints_fraction_value{end+1} = [value1,value2];
                        else
                            error('ERROR - Constraints - Fraction - Minimum fraction must be <= maximum fraction')
                        end
                    else
                        error('ERROR - Constraints - Fraction - Constraint fractions must be >= 0 and <= 1')
                    end
                else
                    error('ERROR - Constraints - Fraction - Constraint values must be a number')
                end
            else
                error('ERROR - Constraints - Fraction - Can only choose one constraint type at a time')
            end
        elseif ~isnan(value1) || ~isnan(value2)
            error('ERROR - Constraints - Fraction - If choosing "Range" constraint, must enter both minimum and maximum fractions');
        end

        % reaction formula
        if added
            value = raw(row,6);
            value = value{1};
            if ~isnan(value)
                if ischar(value)
                    constraints_fraction_reaction{end+1} = value;
                else
                    error('ERROR - Constraints - Fraction - Reaction formula must be valid string')
                end
            else
                error('ERROR - Constraints - Fraction - Must give reaction formula');
            end
        end
        
        % constraints counter
        if added
            constraints_fraction_count = constraints_fraction_count + 1;
        end
    end
end