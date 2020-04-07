
function [constraints_proliferation,constraints_proliferation_type,constraints_proliferation_value] = parse_constraints_proliferation(input_file)
   
    % load data
    [~,~,raw] = xlsread(input_file,'Constraints - Proliferation');

    % proliferation constraint
    row = 3;
    
    % constraint type - equal
    value = raw(row,1);
    value = value{1};
    if ~isnan(value)
        if isfloat(value)           
            if value >= 0
                constraints_proliferation = 1;
                constraints_proliferation_type = 'equal';
                constraints_proliferation_value = value;
            else
                error('ERROR - Constraints - Proliferation - Constraint fraction must be >= 0')
            end
        else
            error('ERROR - Constraints - Proliferation - Constraint fraction must be a number')
        end
    end

    % constraint type - less than
    value = raw(row,2);
    value = value{1};
    if ~isnan(value)
        if ~exist('constraints_proliferation_type','var')
            if isfloat(value)           
                if value >= 0
                    constraints_proliferation = 1;
                    constraints_proliferation_type = 'less';
                    constraints_proliferation_value = value;
                else
                    error('ERROR - Constraints - Proliferation - Constraint fraction must be >= 0')
                end
            else
                error('ERROR - Constraints - Proliferation - Constraint fraction must be a number')
            end
        else
            error('ERROR - Constraints - Proliferation - Can only choose one constraint type')
        end
    end

    % constraint type - greater than
    value = raw(row,3);
    value = value{1};
    if ~isnan(value)
        if ~exist('constraints_proliferation_type','var')
            if isfloat(value)           
                if value >= 0
                    constraints_proliferation = 1;
                    constraints_proliferation_type = 'greater';
                    constraints_proliferation_value = value;
                else
                    error('ERROR - Constraints - Proliferation - Constraint fraction must be >= 0')
                end
            else
                error('ERROR - Constraints - Proliferation - Constraint fraction must be a number')
            end
        else
            error('ERROR - Constraints - Proliferation - Can only choose one constraint type')
        end
    end

    % constraint type - range around - fraction
    value1 = raw(row,4);
    value1 = value1{1};
    value2 = raw(row,5);
    value2 = value2{1};
    if ~isnan(value1) && ~isnan(value2)
        if ~exist('constraints_proliferation_type','var')
            if isfloat(value1) && isfloat(value2)       
                if (value1 >= 0) && (value2 >= 0)
                    constraints_proliferation = 1;
                    constraints_proliferation_type = 'range|fraction';
                    constraints_proliferation_value = [value1,value2];
                else
                    error('ERROR - Constraints - Proliferation - Constraint median and range must both be >= 0')
                end
            else
                error('ERROR - Constraints - Proliferation - Constraint median and range must both be a number')
            end
        else
            error('ERROR - Constraints - Proliferation - Can only choose one constraint type')
        end
    elseif ~isnan(value1) || ~isnan(value2)
        error('ERROR - Constraints - Proliferation - If choosing "Range Around" constraint, must enter both median and range/distance values');
    end

    % constraint type - range around - distance
    value1 = raw(row,6);
    value1 = value1{1};
    value2 = raw(row,7);
    value2 = value2{1};
    if ~isnan(value1) && ~isnan(value2)
        if ~exist('constraints_proliferation_type','var')
            if isfloat(value1) && isfloat(value2)       
                if (value1 >= 0) && (value2 >= 0)
                    constraints_proliferation = 1;
                    constraints_proliferation_type = 'range|distance';
                    constraints_proliferation_value = [value1,value2];
                else
                    error('ERROR - Constraints - Proliferation - Constraint median and range must both be >= 0')
                end
            else
                error('ERROR - Constraints - Proliferation - Constraint median and range must both be a number')
            end
        else
            error('ERROR - Constraints - Proliferation - Can only choose one constraint type')
        end
    elseif ~isnan(value1) || ~isnan(value2)
        error('ERROR - Constraints - Proliferation - If choosing "Range Around" constraint, must enter both median and range/distance values');
    end

    % none
    if ~exist('constraints_proliferation_type','var')
        constraints_proliferation = 0;
        constraints_proliferation_type = NaN;
        constraints_proliferation_value = NaN;
    end
end