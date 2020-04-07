
function [constraints_maxmulti_weight,constraints_maxmulti_reaction,constraints_maxmulti_fraction,constraints_maxmulti_count] = parse_constraints_maxmulti(input_file)
  
    % load data
    [~,~,raw] = xlsread(input_file,'Constraints - Max Multi');

    % initialize max multi constraints
    constraints_maxmulti_weight = {};
    constraints_maxmulti_reaction = {};
    constraints_maxmulti_count = 0;
    constraints_maxmulti_fraction = NaN;

    for row = 2:size(raw,1)

        % weight
        value = raw(row,1);
        value = value{1};
        if ~isnan(value)
            if isfloat(value)
                if value>0

                    % set value
                    constraints_maxmulti_weight{end+1} = value;

                    % reaction
                    value = raw(row,2);
                    value = value{1};
                    if ~isnan(value)
                        if ischar(value)
                            constraints_maxmulti_reaction{end+1} = value;
                            constraints_maxmulti_count = constraints_maxmulti_count + 1;
                        else
                            error('ERROR - Constraints - Max Multi - Reaction formula must be valid string')
                        end
                    else
                        error('ERROR - Constraints - Max Multi - Must give reaction formula');
                    end

                    % fraction of maximum
                    value = raw(1,4);
                    value = value{1};
                    if ~isnan(value)
                        if isfloat(value)
                            if (value >= 0) && (value <= 1)
                                constraints_maxmulti_fraction = value;
                            else
                                error('ERROR - Constraints - Max Multi - Fraction of maximum objective value must either be >= 0 and <= 1')
                            end
                        else
                            error('ERROR - Constraints - Max Multi - Fraction of maximum objective value must either be >= 0 and <= 1')
                        end
                    else
                        error('ERROR - Constraints - Max Multi - Must enter fraction of maximum objective value')
                    end

                else
                    error('ERROR - Constraints - Max Multi - Weight must be > 0')
                end
            elseif isfloat(str2num(value))
                if str2num(value)>0

                    % set value
                    constraints_maxmulti_weight{end+1} = str2num(value);

                    % reaction
                    value = raw(row,2);
                    value = value{1};
                    if ~isnan(value)
                        if ischar(value)
                            constraints_maxmulti_reaction{end+1} = value;
                            constraints_maxmulti_count = constraints_maxmulti_count + 1;
                        else
                            error('ERROR - Constraints - Max Multi - Reaction formula must be valid string')
                        end
                    else
                        error('ERROR - Constraints - Max Multi - Must give reaction formula');
                    end

                    % fraction of maximum
                    value = raw(1,4);
                    value = value{1};
                    if ~isnan(value)
                        if isfloat(value)
                            if (value >= 0) && (value <= 1)
                                constraints_maxmulti_fraction = value;
                            else
                                error('ERROR - Constraints - Max Multi - Fraction of maximum objective value must either be >= 0 and <= 1')
                            end
                        else
                            error('ERROR - Constraints - Max Multi - Fraction of maximum objective value must either be >= 0 and <= 1')
                        end
                    else
                        error('ERROR - Constraints - Max Multi - Must enter fraction of maximum objective value')
                    end

                else
                    error('ERROR - Constraints - Max Multi - Weight must be > 0')
                end
            else
                error('ERROR - Constraints - Max Multi - Weight must be a number')
            end
        end
    end
end