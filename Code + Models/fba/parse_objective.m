
function [objective_function_direction,objective_function_weight,objective_function_reaction] = parse_objective(input_file)
  
% load data
[~,~,raw] = xlsread(input_file,'Objective Function');

% initialize objective functions
objective_function_direction = {};
objective_function_weight = {};
objective_function_reaction = {};

for row = 2:size(raw,1)

    % objective function
    value = raw(row,1);
    value = value{1};
    if ~isnan(value)
        if ischar(value)
            if strcmpi(value,'X')

                % max or min
                value = raw(row,2);
                value = value{1};
                if ~isnan(value)
                    if ischar(value)
                        if any(strcmp({'MAX','MIN'},value))
                            objective_function_direction{end+1} = value;
                        else
                            error('ERROR - Objective Function - Must state MAX or MIN')
                        end
                    else
                        error('ERROR - Objective Function - Must state MAX or MIN')
                    end
                else
                    error('ERROR - Objective Function - Must state MAX or MIN');
                end

                % weight
                value = raw(row,3);
                value = value{1};
                if ~isnan(value)
                    if isfloat(value)
                        if value>0
                            objective_function_weight{end+1} = value;
                        else
                            error('ERROR - Objective Function - Weight must be > 0')
                        end
                    elseif isfloat(str2num(value))
                        if str2num(value)>0
                            objective_function_weight{end+1} = str2num(value);
                        else
                            error('ERROR - Objective Function - Weight must be > 0')
                        end
                    else
                        error('ERROR - Objective Function - Weight must be a number')
                    end
                else
                    error('ERROR - Objective Function - Weight must be a number');
                end

                % reaction formula
                value = raw(row,4);
                value = value{1};
                if ~isnan(value)
                    if ischar(value)
                        objective_function_reaction{end+1} = value;
                    else
                        error('ERROR - Objective Function - Reaction formula must be valid string')
                    end
                else
                    error('ERROR - Objective Function - Must give reaction formula');
                end

           else
                error('ERROR - Objective Function - Entry in SELECT column must either be X or blank')
            end
        else
            error('ERROR - Objective Function - Entry in SELECT column must either be X or blank')
        end
    end
end