
function [use_mutation] = parse_mutations(input_file)

    % load data
    [~,~,raw] = xlsread(input_file,'Mutations');
    
    % use mutation data
    row = 3;
    value = raw(row,1);
    value = value{1};
    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                use_mutation = 1;
            else
                error('ERROR - Mutations - Value must either be X or blank')
            end
        else
            error('ERROR - Mutations - Value must either be X or blank')
        end
    end
    
    row = 4;
    value = raw(row,1);
    value = value{1};
    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                if ~exist('use_mutation','var')
                    use_mutation = 2;
                else
                    error('ERROR - Mutations - Cannot select more than one option') 
                end
            else
                error('ERROR - Mutations - Value must either be X or blank')
            end
        else
            error('ERROR - Mutations - Value must either be X or blank')
        end
    end
    
    row = 5;
    value = raw(row,1);
    value = value{1};
    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')
                if ~exist('use_mutation','var')
                    use_mutation = 0;
                else
                    error('ERROR - Mutations - Cannot select more than one option') 
                end
            else
                error('ERROR - Mutations - Value must either be X or blank')
            end
        else
            error('ERROR - Mutations - Value must either be X or blank')
        end
    end
    
    if ~exist('use_mutation','var')
        error('ERROR - Mutations - Must select an option') 
    end
end